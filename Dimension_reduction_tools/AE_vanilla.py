from tensorflow.keras import Input, layers, losses, metrics
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import optimizers

import numpy as np
import numpy.linalg as la 
import scipy.io
import matplotlib.pyplot as plt

# load training and test sets
mat = scipy.io.loadmat('X')
x_train = np.transpose(mat['X'])
mat_1 = scipy.io.loadmat('Y')
x_test = np.transpose(mat_1['Y'])

print('shape is',x_train.shape)


# degree of freedoms for CVs
encoding_dim = 2  


# set Input
input_img = Input(shape=(x_train.shape[1],))

# layers between input and bottleneck
encoded = layers.Dense(1000, activation='tanh')(input_img)
encoded_2 = layers.Dense(encoding_dim, activation='tanh')(encoded)

# set encoder
encoder = Model(input_img, encoded_2)

# set input for decoder
decoder_input = Input(shape=(encoding_dim,))

# Layers between bottleneck and output
decoded = layers.Dense(1000, activation='tanh')(decoder_input)
decoded_2 = layers.Dense(x_train.shape[1])(decoded)

# set decoder
decoder = Model(decoder_input,decoded_2)

# 
auto_encoder_input =  Input(shape=(x_train.shape[1],))
auto_encoder_encoder_out = encoder(auto_encoder_input)
auto_encoder_decoder_out = decoder(auto_encoder_encoder_out)
autoencoder = Model(auto_encoder_input, auto_encoder_decoder_out)

callback = EarlyStopping(monitor='val_loss', patience=50)


opt = optimizers.SGD(learning_rate=0.01)
autoencoder.compile(optimizer=opt, loss='mse')


history = autoencoder.fit(x_train, x_train,
                epochs=2000,
                batch_size=20,
                shuffle=True,
                callbacks=[callback],
                verbose=1,
                validation_data=(x_test, x_test))
# loss = 1/(M*d)* sum of squares of errors, where M is the size of training set and d is the dimension of input vector
# validation loss is defined similarly over the test set, instead of M, m as the size of test set

# Encoder and decoder
encoded_imgs = encoder.predict(x_train)
decoded_imgs = decoder.predict(encoded_imgs)

CV = encoder.predict(x_train)

file_name = 'CV.mat'
scipy.io.savemat(file_name, {'CV': CV})


plt.plot(history.history["loss"], label="Training Loss")
plt.plot(history.history["val_loss"], label="Validation Loss")
plt.legend()
plt.savefig('re01.png')
print('loss',x_train.shape[1]*history.history["loss"][-1])
print('Val_loss',x_train.shape[1]*history.history["val_loss"][-1])
print('relative_error',x_train.shape[0]*x_train.shape[1]*history.history["loss"][-1]/la.norm(x_train,'fro'))
