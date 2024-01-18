import numpy as np
from matplotlib import pyplot # For displaying results
from typing import List

def show_image(X: np.ndarray, y: str) -> None:
    """
    Shows an image using matplotlib

    Arguments:
        - X (np.ndarray) - image
        - y (str) - image label
    """

    pyplot.imshow(X, cmap=pyplot.get_cmap('gray'))
    pyplot.title(f'Number {y}')
    pyplot.show()

def to_vectors(y: List[int]) -> np.ndarray:
    """
    Convert a list of labels to the vector format

    Example:
    to_vectors([1, 2]) = [[0,1,0,...,0], [0,0,1,0,...,0]]
    
    Arguments: 
        y (List[int]) - a list of labels
    """

    return np.array([[1.0 if i == label else 0.0 for i in range(10)] for label in y])
