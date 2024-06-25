#! usr/bin/env python

import numpy as np

class DecisionTreeClassifier:
    def __init__(self, max_depth=0):
        self.max_depth = max_depth # Max depth of tree.
        self.__fit_data = None # To be set by fit(X, y).

    def __entropy(self, y):
        # Return 0 entropy for an empty array.
        if len(y) == 0:
            return 0
        
        # Sum probability * log2(probability) for each label.
        unique, counts = np.unique(y, return_counts=True)
        probabilities = counts / len(y)
        bits = np.sum(-probabilities * np.log2(probabilities))

        return bits
    
    def __info_gain(self, X, y, threshold, feature_index):
        # Determine splits according to threshold value.
        feature = X[:,feature_index]

        left_indices = np.argwhere(feature < threshold).flatten()
        right_indices = np.argwhere(feature >= threshold).flatten()
        
        # Split y.
        left_child = y[left_indices]
        right_child = y[right_indices]

        # Calculate probabilities for child groups.
        prob_left = left_child.shape[0] / y.shape[0]
        prob_right = right_child.shape[0] / y.shape[0]

        info_gain = self.__entropy(y) - (prob_left * self.__entropy(left_child) + prob_right * self.__entropy(right_child))
        return info_gain
    
    def __find_best_split(self, X, y):
        best_feature = (0, None, None)

        # Iterate over features.
        for feature_index in range(X.shape[1]):
            feature = X[:,feature_index]
            feature_values = np.unique(feature)
            
            # Iterate over thresholds.
            for threshold in feature_values:
                info = self.__info_gain(X, y, threshold, feature_index)

                # Save feature and threshold with highest information gain.
                if info > best_feature[0]:
                    best_feature = (info, feature_index, threshold)

        return best_feature[1], best_feature[2]

    def __split_dataset(self, X, y, feature_index, threshold):
        # Determine splits according to threshold value.
        feature = X[:,feature_index]
        left_indices = np.argwhere(feature < threshold).flatten()
        right_indices = np.argwhere(feature >= threshold).flatten()
        
        # Split X and y.
        left_X, left_y = X[left_indices], y[left_indices]
        right_X, right_y = X[right_indices], y[right_indices]

        return left_X, left_y, right_X, right_y

    class __DecisionTreeNode:
        def __init__(self, feature_index=None, threshold=None, left_child=None, right_child=None, label=None):
            self.feature_index = feature_index  # Index of the feature for splitting
            self.threshold = threshold  # Threshold for splitting
            self.left_child = left_child  # Left child node
            self.right_child = right_child  # Right child node
            self.label = label  # Label for leaf nodes

    def __build_tree(self, X, y, max_depth):
        # Initialize root node.
        node = self.__DecisionTreeNode()
        
        # Return leaf node with label at end of recursion.
        if max_depth == 0 or len(np.unique(y)) == 1:
            # If maximum depth reached or all labels are the same
            node.label = np.argmax(np.bincount(y)) # Most common label
            return node
        
        # If not leaf node, split into left and right children and continue recursion.
        else:
            split_feature, split_threshold = self.__find_best_split(X, y)
            left_X, left_y, right_X, right_y = self.__split_dataset(X, y, split_feature, split_threshold)
            
            node.feature_index = split_feature
            node.threshold = split_threshold
            
            # Recursive call to build left and right child nodes
            node.left_child = self.__build_tree(left_X, left_y, max_depth - 1)
            node.right_child = self.__build_tree(right_X, right_y, max_depth - 1)
    
        return node

    def __predict_tree(self, X, tree):
        y = []
        
        # For each point, go through beginning of tree until leaf node is reached
        for point in X:
            node = tree
            while node.label is None:
                if point[node.feature_index] < node.threshold:
                    node = node.left_child
                else:
                    node = node.right_child

            # Use leaf node label for predicted y.
            y.append(node.label)
        return np.array(y)

    def fit(self, X, y):
        self.__fit_data = self.__build_tree(X, y, self.max_depth)
        return
    
    def predict(self, X):
        return self.__predict_tree(X, self.__fit_data)

def main():
    from sklearn import datasets
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error

    # Load dataset.
    iris = datasets.load_iris()

    # Split into training and test sets.
    X = iris.data
    y = iris.target
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1)

    # Initialize decision tree.
    tree = DecisionTreeClassifier(max_depth = 2)
    tree.fit(X_train, y_train)

    # Compare predicted and actual labels.
    mse = mean_squared_error(y_test, tree.predict(X_test))
    print(f'MSE: {mse}')

    return mse

if __name__ == "__main__":
    main()