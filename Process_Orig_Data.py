import numpy as np
from scipy.interpolate import interp1d
from Savitzky_Golay_filter import savitzky_golay
from Lowess_filter import lowess

class process_orig_data:
    def __init__(self, input_file):
        self.input_file = input_file
        self.prefix = self.input_file.split('.')[1]
        self.scaling_factor = 1
        self.get_scaling_factor()
        self.useful_data = np.array([])
        self.extract_useful_data()
        self.plateau_data = np.array([])
        self.extract_plateau()
        self.L = 0
        self.get_L()
        self.smooth_data_lowess = np.array([])
        self.get_smooth_data_lowess()

    def get_scaling_factor(self):
        fin = open(self.input_file, 'r')
        for line in fin:
            if line.startswith("Dacto[A]xy"):
                self.scaling_factor = float(line.split('=')[1])
                break
        fin.close()

    def extract_useful_data(self):
    ##-----extract the matrix data points at the end of file-----#
        with open(self.input_file) as textFile:
            gaga = textFile.read().splitlines()
        ##----- maybe different for different measurements files!!
        data_in_the_end = np.array([[float(i) for i in line.split()] for line in gaga[575:] ])
        self.useful_data = data_in_the_end[:, [0, 1, 2, 4]]

    def extract_plateau(self):
    ##-----look for the endpoints of plateau-----##
        cutoff = 1.5
        st_point = np.argmax(self.useful_data[:,3])
        while abs(self.useful_data[st_point-1, 3] - self.useful_data[st_point, 3]) < cutoff:
            st_point -= 1
        left_end = st_point + 1
        st_point = int(round(self.useful_data.shape[0]/2))
        while abs(self.useful_data[st_point+1, 3] - self.useful_data[st_point, 3]) < cutoff:
            st_point += 1
        right_end = st_point - 1
        self.plateau_data = self.useful_data[left_end:right_end+1,:]

    def get_L(self):
    ##-----calculate the parameter L-----##
        X_Y = self.plateau_data[:,1:3]
        X_Y = X_Y - np.ones_like(X_Y) * X_Y[0]
        self.L = self.scaling_factor * np.sqrt( np.square(X_Y[:,0]) + np.square(X_Y[:,1]) ) *3

    def get_smooth_data_lowess(self):
    ##-----perform Lowess filter, output is a n x 2 matrix-----##
        smooth_data_lowess = lowess(self.L, self.plateau_data[:,3], 0.3, 3)
        self.smooth_data_lowess = np.concatenate((self.L[:, np.newaxis], smooth_data_lowess[:, np.newaxis]), axis = 1)
