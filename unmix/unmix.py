import tensorflow as tf
import scipy.io as cpio
import numpy as np

dt = np.float32

# configuration
iterations = 5000

# load data to unmix
mat = cpio.loadmat('temp.mat')
Y = mat['Y'].astype(dt)
waveform = np.squeeze(mat['waveform'].astype(dt))

# get dimensions
components = Y.shape[0]
timesteps = Y.shape[1]
waveform_len = len(waveform)

# M0?
if 'M0' in mat:
    M0 = mat['M0'].astype(dt)
    assert(M0.shape == (components, components))
else:
    M0 = np.eye(components, components, dtype=dt) + 0.00001
    print('Seeding random M_0')
if 'X0' in mat:
    X0 = mat['X0'].astype(dt)
    assert(X0.shape == (components, timesteps + waveform_len - 1))
else:
    X0 = np.random.rand(components, timesteps + waveform_len - 1).astype(dt)
    print('Seeding random X_0')

# pass Y to tensorflow
tf_Y = tf.constant(Y)
tf_waveform = tf.constant(waveform[::-1])

# unmixing components
tf_M = tf.Variable(M0, constraint=lambda v: tf.clip_by_value(v, 0, np.inf))
tf_X = tf.Variable(X0, constraint=lambda v: tf.clip_by_value(v, 0, np.inf))

# convolve
tf_X4D = tf.expand_dims(tf.expand_dims(tf_X, 1), -1)
tf_filter = tf.reshape(tf_waveform, (1, -1, 1, 1))
tf_Xw = tf.squeeze(tf.nn.conv2d(tf_X4D, tf_filter, strides=[1, 1, 1, 1], padding="SAME"))

tf_MX = tf.matmul(tf_M, tf_Xw)

tf_sparse_M_row = (np.sqrt(components) - tf.norm(tf_M, ord=1, axis=0) / (1e-5 + tf.norm(tf_M, ord=2, axis=0))) / \
              (np.sqrt(components) - 1)

tf_sparse_X_row = (np.sqrt(components) - tf.norm(tf_X, ord=1, axis=0) / (1e-5 + tf.norm(tf_X, ord=2, axis=0))) / \
              (np.sqrt(components) - 1)

tf_sparse_X_col = (np.sqrt(timesteps) - tf.norm(tf_X, ord=1, axis=1) / (1e-5 + tf.norm(tf_X, ord=2, axis=1))) / \
              (np.sqrt(timesteps) - 1)

# cost
tf_err = tf_Y - tf_MX[:, waveform_len-1:]
# tf_err = tf.divide(tf_Y - tf_MX[:, waveform_len-1:], tf_MX[:, waveform_len-1:])
tf_cost = tf.reduce_sum(tf.pow(tf_err, 2))

# steps
opt = tf.train.AdamOptimizer()
step_train = opt.minimize(tf_cost)

step_init = tf.global_variables_initializer()

# track costs
cost = np.nan
costs = np.zeros((iterations, ))
stop_condition = 1  # max iterations

with tf.Session() as sess:
    sess.run(step_init)
    for i in range(iterations):
        # update X and get cost
        _, cost = sess.run((step_train, tf_cost))

        # store cost
        costs[i] = cost

        # log every so often
        if 0 == (i % 100):
            print('Iteration #{}: {}'.format(i + 1, cost))

    M = sess.run(tf_M)
    Xw = sess.run(tf_Xw)
    X = sess.run(tf_X)

    # trim
    Xw = Xw[:, waveform_len-1:]
    X = X[:, waveform_len-1:]

# write output
cpio.savemat('temp.mat', dict([('M', M), ('Xw', Xw), ('X', X), ('cost', cost), ('costs', costs),
                               ('stop_condition', stop_condition)]))
