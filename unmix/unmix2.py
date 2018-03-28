import tensorflow as tf
import scipy.io as cpio
import numpy as np
import math

dt = np.float32

stop_threshold = 1e-5
stop_without_improvement = 500

# load data to unmix
mat = cpio.loadmat('temp.mat')
Y = mat['Y'].astype(dt)
waveform = np.squeeze(mat['waveform'].astype(dt))

# get dimensions
components = Y.shape[0]
timesteps = Y.shape[1]
waveform_len = len(waveform)

# pass Y to tensorflow
tf_Y = tf.constant(Y)
tf_waveform = tf.constant(waveform[::-1])

# unmixing components
tf_M = tf.Variable(np.random.rand(components, components).astype(dt))
tf_X = tf.Variable(np.random.rand(components, timesteps + waveform_len - 1).astype(dt))

# convolve
tf_X4D = tf.expand_dims(tf.expand_dims(tf_X, 1), -1)
tf_filter = tf.reshape(tf_waveform, (1, -1, 1, 1))
tf_Xw = tf.squeeze(tf.nn.conv2d(tf_X4D, tf_filter, strides=[1, 1, 1, 1], padding="SAME"))

tf_MX = tf.matmul(tf_M, tf_Xw)

# error normalized
tf_err = tf_Y - tf_MX[:, waveform_len-1:]
#tf_err = tf.divide(tf_Y - tf_MX[:, waveform_len-1:], tf_MX[:, waveform_len-1:])

# cost
sps = 20
freq = 0.4
p = freq / sps
prob_M = tf.distributions.Exponential(0.00638246)
prob_X = tf.distributions.Normal(float(timesteps) * p, math.sqrt(float(timesteps) * p * (1 - p)))
prob_err = tf.distributions.Normal(0.0, 0.1)

tf_prob = tf.reduce_sum(prob_M.log_prob(tf_M)) + tf.reduce_sum(prob_X.log_prob(tf_X)) + tf.reduce_sum(prob_err.log_prob(tf_err))
tf_cost = tf.reduce_sum(tf.pow(tf_Y - tf_MX[:, waveform_len-1:], 2))

# configuration
iterations = 100000
opt = tf.train.GradientDescentOptimizer(0.001)
gradients = opt.compute_gradients(0 - tf_prob)
clipped_gradients = [(tf.clip_by_value(grad, -np.inf, var), var) for grad, var in gradients]
step_train = opt.apply_gradients(clipped_gradients)
step_init = tf.global_variables_initializer()

# force nonnegativity via clip
clip_M = tf_M.assign(tf.multiply(tf_M, tf.cast(tf_M >= 0, dtype=dt)))
clip_X = tf_X.assign(tf.multiply(tf_X, tf.cast(tf_X >= 0, dtype=dt)))
step_clip = tf.group(clip_M, clip_X)

# minimum
stop_threshold *= components * components + components * timesteps
mn = np.inf
no_improvement = 0
stop_condition = 1  # max iterations
cost = None
costs = np.zeros((iterations, ))

with tf.Session() as sess:
    sess.run(step_init)
    for i in range(iterations):
        _, cost, prob = sess.run((step_train, tf_cost, tf_prob))
        # sess.run(step_clip)
        costs[i] = cost

        if np.isnan(cost):
            raise RuntimeError('Cost became NaN')

        if cost < mn:
            mn = cost
            no_improvement = 0
        elif cost < stop_threshold:
            stop_condition = 2  # below threshold
            costs = costs[0:iterations+1]
            break
        else:
            no_improvement += 1
            if no_improvement > stop_without_improvement:
                stop_condition = 3  # no improvement
                costs = costs[0:iterations + 1]
                break

        if 0 == (i % 100):
            print('Iteration #{}\t{:.2f}\t{:.2f}'.format(i + 1, cost, prob))

    M = sess.run(tf_M)
    Xw = sess.run(tf_Xw)
    X = sess.run(tf_X)

    # trim
    Xw = Xw[:, waveform_len-1:]
    X = X[:, waveform_len-1:]

# write output
cpio.savemat('temp.mat', dict([('M', M), ('Xw', Xw), ('X', X), ('cost', cost), ('costs', costs), ('stop_condition', stop_condition)]))
