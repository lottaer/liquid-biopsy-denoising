#     Estimation of purity and subclonal ratio from liquid biosy sequencing. 
#     Inspired by Lakatos et. al., <https://github.com/elakatos/liquidCNA>.
#     Copyright (C) 2024  Lotta Eriksson & Linnea Hallin
#                    lottaer@chalmers.se   hallinl@chalmers.se
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.

import platform
import os
import keras
import pandas as pd
import numpy as np
import tensorflow as tf
import keras.callbacks as callbacks
import keras.layers as layers
from keras.models import Model
from keras.callbacks import History

# Sets the learning rate of the Adam oprimizers and returns it
def get_adam_optimizer(learning_rate = 0.0001):
  return tf.keras.optimizers.Adam(learning_rate = learning_rate)
 
# Print model parameters 
def print_model_summary(model):
  model.summary()
  return None 
  
# Load a saved model from disk
def load_model(path = "denoising_model"):
  return keras.models.load_model(path, compile = False)
  
# Reshape R-objects and convert to numpy array
def to_reshaped_array(robject):
  robject = np.array(robject)
  object_shape = robject.shape
  return robject.reshape((object_shape[0], object_shape[1], 1))

# Custom loss function (MSE + total variation) 
@tf.function
def total_variation_mse(y_true, y_pred, penalty = 0.1):
    mse_loss = tf.reduce_mean((y_pred - y_true) ** 2)
    diff = y_pred[:, 1:, :] - y_pred[:, :-1, :]
    total_var_loss = tf.reduce_mean(tf.abs(diff))
    return mse_loss + penalty * total_var_loss

# Define the model
# segments: the dround truth (columns being genomic bins)
# noisy: noisy version of segments
# penalty: penalty parameters to use in loss function
# n_epochs: number of epochs the model is trained
# feature_size: number of feature maps in the first layer
# n_layers: number of hidden layers
# kernel_size: size of the kernel
# pooling_method: pooling method to use, if not 'max' AveragePooling is used
# save_path: location to save the trained model
def run_conv1d_denoiser(segments, noisy, penalty = 0.025, n_epochs = 10,
                        feature_size = 128, n_layers = 2, kernel_size = 7,
                        pooling_method = 'max', save_path = None):
  # convert R data into numpy array
  original_data = to_reshaped_array(segments)
  corrupted_data = to_reshaped_array(noisy)
  input_shape = (corrupted_data.shape[1], 1)
  n_layers = int(n_layers)
  feature_size = int(feature_size)
  kernel_size = int(kernel_size)
  def tvmse(y_true, y_pred): 
    return total_variation_mse(y_true, y_pred, penalty)
  # define the model
  model = keras.Sequential()
  model.add(layers.Conv1D(feature_size, kernel_size = kernel_size, 
                          activation = "relu", padding = "same"))
  # encoder part
  for i in range(n_layers):
    feature_size /= 2 
    feature_size = int(feature_size)
    model.add(layers.Conv1D(feature_size, kernel_size = kernel_size,
                          activation = "relu", padding = "same"))
    if pooling_method == 'max':
      model.add(layers.MaxPooling1D(2))
    else: 
      model.add(layers.AveragePooling1D(2))
  # decoder part
  for i in range(n_layers):
    feature_size *= 2 
    feature_size = int(feature_size)
    model.add(layers.Conv1DTranspose(feature_size, kernel_size = kernel_size,
                          activation = "relu", padding = "same", 
                          strides = 1))
    model.add(layers.UpSampling1D(2))
  # add the output layer
  model.add(layers.Conv1DTranspose(1, activation = "relu",
                                   padding = "same", kernel_size = 1,
                                   strides = 1))
  opt = get_adam_optimizer()
  model.compile(optimizer = opt, loss = tvmse, run_eagerly = False)
  stop_early = tf.keras.callbacks.EarlyStopping(monitor = "val_loss", 
                  patience = 2)
  history = model.fit(corrupted_data, original_data, shuffle = True, 
                      epochs = n_epochs, validation_split = 0.2,
                      callbacks = [stop_early])
  # save model to not have to train it again
  if save_path is not None:
    model.save(save_path)
  return {"model": model, "history": history.history}

# Predicts a new observation
def predict_observations(model, new_data):
  new_data = np.array(new_data) # convert Robject to numpy array
  if new_data.ndim == 1: 
    new_data = new_data[None, ...]
    # print(new_data.shape)
  return model.predict(new_data, verbose = 0)
