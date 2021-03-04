import numpy as np
import matplotlib.pyplot as plt

def plot_csv(file='tagawa.csv'):
    Ha, maxU, gradU = np.loadtxt(file, unpack=True)
    maxU = np.array(maxU)
    gradU = np.array(gradU)
    Ha = ['Ha='+str(int(ha)) for ha in Ha]
    true_maxU = np.array([16.43, 4.073, 1.123])
    true_gradU = np.array([50.00, 10.0, 2.5])
    rel_error_grad = abs((gradU-true_gradU)/true_gradU)*100
    rel_error_maxU = abs((maxU-true_maxU)/true_maxU)*100

    fig_gradU, ax_gradU = plt.subplots(figsize=(12,12))
    X = np.arange(3) # for ax.bar() purposes
    ax_gradU.bar(X + 0.0, gradU, color = 'b', width = 0.25, label='Analytical Gradient')
    ax_gradU.bar(X + 0.25, true_gradU, color = 'g', width = 0.25, label='Q2D Gradient')
    ax_gradU.bar(X + 0.5, rel_error_grad, color = 'r', width = 0.25, label='relative error')
    ax_gradU.set_xticks(X)
    ax_gradU.set_xticklabels(Ha)
    ax_gradU.set_title('Maximum Velocity')

    fig_maxU, ax_maxU = plt.subplots(figsize=(12,12))
    X = np.arange(3) # for ax.bar() purposes
    ax_maxU.bar(X + 0.0, maxU, color = 'b', width = 0.25, label='tagawa max U')
    ax_maxU.bar(X + 0.25, true_maxU, color = 'g', width = 0.25, label='Q2D max U')
    ax_maxU.bar(X + 0.5, rel_error_grad, color = 'r', width = 0.25, label='relative error')
    ax_maxU.set_xticks(X)
    ax_maxU.set_xticklabels(Ha)
    ax_maxU.set_title('Maximum Velocity')

    ax_gradU.grid()
    ax_gradU.legend(loc='upper right')
    ax_maxU.grid()
    ax_maxU.legend(loc='upper right')
    #ax.bar(X + 0.25, data[1], color = 'g', width = 0.25)
    #ax.bar(X + 0.50, data[2], color = 'r', width = 0.25)

    plt.show()

plot_csv()
