#! usr/bin/env python

import numpy as np

class LSTM:
    def __init__(self, input_size, hidden_size):
        self.input_size = input_size
        self.hidden_size = hidden_size

        # Initialize weights for input.
        self.Ui = np.random.randn(hidden_size, input_size)
        self.Uf = np.random.randn(hidden_size, input_size)
        self.Uc = np.random.randn(hidden_size, input_size)
        self.Uo = np.random.randn(hidden_size, input_size)

        # Initialize weights for previous hidden layer.
        self.Wi = np.random.randn(hidden_size, hidden_size)
        self.Wf = np.random.randn(hidden_size, hidden_size)
        self.Wc = np.random.randn(hidden_size, hidden_size)
        self.Wo = np.random.randn(hidden_size, hidden_size)

        # Initialize weights for previous memory cell.
        self.Pi = np.random.randn(hidden_size, hidden_size)
        self.Pf = np.random.randn(hidden_size, hidden_size)
        self.Po = np.random.randn(hidden_size, hidden_size)

    def sigmoid(self, x):
        return 1 / (1 + np.exp(-x))

    def tanh(self, x):
        return np.tanh(x)

    def forward(self, xt, ht_1, ct_1):
        # Input gate.
        it = self.sigmoid(self.Ui @ xt + self.Wi @ ht_1 + self.Pi @ ct_1)
        
        # Forget gate.
        ft = self.sigmoid(self.Uf @ xt + self.Wf @ ht_1 + self.Pf @ ct_1)
        
        # Update gate.
        ctilde_t = self.tanh(self.Uc @ xt + self.Wc @ ht_1)
        ct = ft * ct_1 + it * ctilde_t
        
        # Output gate.
        ot = self.sigmoid(self.Uo @ xt + self.Wo @ ht_1 + self.Po @ ct)
        ht = ot * self.tanh(ct)
        
        return ht, ct

def forward_pass(input_size, hidden_size, number_passes):
    lstm = LSTM(input_size, hidden_size)
    x = np.random.randn(3)
    h = np.random.randn(4)
    c = np.random.randn(4)
    for _ in range(number_passes):
        h, c = lstm.forward(x, h, c)
    print("Output hidden state:", h)
    print("Output memory cell:", c)
    return h,c

if __name__ == "__main__":
    print(forward_pass(3, 4, 10))