from keras.models import Sequential
from keras.layers import Dense
import kerastuner as kt

import pandas as pd
import glob
import sys
import numpy as np
from matplotlib import pyplot

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

HIDDEN_SIZE = 800

STRANDS    = ['+', '-']
CHROMOS    = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
SV_CLASSES = ['TRA', 'DEL', 'DUP', 'h2hINV', 't2tINV']
SV_METHODS = ['SNOWMAN_dRANGER', 'SNOWMAN_BRASS', 'SNOWMAN_dRANGER_DELLY', 'BRASS_dRANGER_DELLY',
              'SNOWMAN_DELLY', 'dRANGER_DELLY', 'SNOWMAN_BRASS' 'BRASS_DELLY',
              'SNOWMAN_BRASS_DELLY', 'BRASS_dRANGER', 'SNOWMAN_BRASS_dRANGER_DELLY',
              'SNOWMAN_BRASS_dRANGER', 'BRASS_DELLY']

def encode(df):
    df['chrom1']    = df['chrom1'].map(lambda a: CHROMOS.index(a))
    df['chrom2']    = df['chrom2'].map(lambda a: CHROMOS.index(a))
    df['strand1']   = df['strand1'].map(lambda a: STRANDS.index(a))
    df['strand2']   = df['strand2'].map(lambda a: STRANDS.index(a))
    df['svclass']   = df['svclass'].map(lambda a: SV_CLASSES.index(a)) # +1?
    df['svmethod']  = df['svmethod'].map(lambda a: SV_METHODS.index(a))
    return df

def load_data(save_path='sv_data.npy', load_path=None):
    if load_path != None:
        print('Loading data from files: {}'.format(load_path))
        return np.load(load_path + '_X.npy'), np.load(load_path + '_y.npy')

    df_target = pd.read_csv('csv/dataframe.csv', sep=',')

    print('Calculating filepaths...')
    path = r'dataset/sv'
    all_files = glob.glob(path + "/*.bedpe")
    li = []
    li_target = []
    lengths = []

    print('Preparing data...')
    for i, filename in enumerate(all_files):
        df = pd.read_csv(filename, sep='\t', dtype={'chrom1': str, 'chrom2': str})
        df.drop(['sv_id'], axis=1, inplace=True)
        df = encode(df)
        lengths.append(df.shape[0])
        df.index = df.index + 1
        df_out = df.stack()
        df_out.index = df_out.index.map('{0[1]}_{0[0]}'.format)
        df_out = df_out.to_frame().T
        cols = df_out.shape[1]
        fn = filename[11:]
        df_out['genome'] = fn[:-44]
        df_out = df_out.merge(df_target, how='left', on='genome')
        df_out.drop(['genome'], axis=1, inplace=True)
        [unpacked] = df_out.to_numpy(dtype=np.int)
        li_target.append(unpacked[cols:])
        li.append(unpacked[:cols])

        sys.stdout.write("\r%d%%" % (float(i)/float(len(all_files)) * 100.0))
        sys.stdout.flush()

    print('Padding data...')
    
    arr = np.zeros([len(li),len(max(li, key=lambda x: len(x)))])
    for i, j in enumerate(li):
        arr[i][0:len(j)] = j

    np.save(save_path + '_X', arr)
    np.save(save_path + '_y', li_target)
    return arr, np.asarray(li_target)

def create_model(X_train, y_train):
    model = Sequential()
    model.add(Dense(HIDDEN_SIZE, input_dim=X_train.shape[1], activation='relu'))
    model.add(Dense(HIDDEN_SIZE, activation='relu'))
    model.add(Dense(y_train.shape[1], activation='sigmoid'))
    model.compile(loss='mse', optimizer='adam', metrics=['binary_accuracy'])
    return model

def main():
    # Load from file or re-generate
    X, y = load_data(load_path='sv_data.npy')
    print('-- Data size --\nX: {} \ny: {}'.format(X.shape, y.shape))

    # Split test and train set
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

    # Scale data
    scaler = MinMaxScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Create the model
    # model = create_model(X_train, y_train)

    # # Train
    # history = model.fit(X_train, y_train, epochs=80, batch_size=10, validation_data=(X_test, y_test), shuffle=True)


    # _, train_accuracy = model.evaluate(X_train, y_train)
    # _, test_accuracy = model.evaluate(X_test, y_test)

    # print('Accuracy (training): %.2f' % (train_accuracy * 100))
    # print('Accuracy (testing): %.2f' % (test_accuracy * 100))

    # # plot loss during training
    # pyplot.subplot(211)
    # pyplot.title('Loss')
    # pyplot.xlabel('epoch')
    # pyplot.plot(history.history['loss'], label='train')
    # pyplot.plot(history.history['val_loss'], label='test')
    # pyplot.legend()

    # # plot accuracy during training
    # pyplot.subplot(212)
    # pyplot.title('Accuracy')
    # pyplot.xlabel('epoch')
    # pyplot.plot(history.history['binary_accuracy'], label='train')
    # pyplot.plot(history.history['val_binary_accuracy'], label='test')
    # pyplot.legend()

    # pyplot.tight_layout()
    # pyplot.show()

    model = create_model(X_train, y_train)


    tuner = kt.Hyperband(model,
                     objective='val_binary_accuracy',
                     max_epochs=10,
                     factor=3,
                     directory='./',
                     project_name='intro_to_kt')


if __name__ == "__main__":
    main()
