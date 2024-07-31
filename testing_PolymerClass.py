import numpy as np
from PolymerClass import Polymer  # Import the Polymer class


def test_unit_basis_3pts():
    """
    START with two vectors in the yz plane pointing in the y direction
    END with unit vectors (k: y[meta CCC arrow], w: x[ortho CCC], v: z [para CCC])
    :return: NULL
    """
    polymer = Polymer(10, [-20, 20, -20, 20, -20, 20])

    # input arrow > along xz
    c1 = np.array([0, -0.5, 0.5])
    c2 = np.array([0.0, 0.0, 0.0])
    c3 = np.array([0.0, -0.5, -0.5])
    left_center_right = [c1, c2, c3]

    # Expected output- unit basis
    k_expected = np.array([0.0, 1.0, 0.0])
    w_expected = np.array([1.0, 0.0, 0.0])
    v_expected = np.array([0.0, 0.0, 1.0])

    # Run the method
    k, w, v = polymer._unit_basis_3pts(left_center_right)

    # Check the results
    assert np.allclose(k, k_expected), f"Expected k: {k_expected}, but got: {k}"
    assert np.allclose(w, w_expected), f"Expected w: {w_expected}, but got: {w}"
    assert np.allclose(v, v_expected), f"Expected v: {v_expected}, but got: {v}"

    print("All tests passed.")

def test_rodrigues_rotation():
    """
    START x axis, ROTATE 90 with z axis
    END with y axis
    :return: NULL
    """
    # input arrow > along xz
    x = np.array([1, 0, 0])
    z = np.array([0, 0, 1])
    y = np.array([0, 1, 0])

    # Run the method
    y_expected = Polymer._rodrigues_rotation(90, x, z)

    # Check the results
    assert np.allclose(y, y_expected), f"Expected k: {y_expected}, but got: {y}"

    print("All tests passed.")


# Run the test
# test_unit_basis_3pts()
test_rodrigues_rotation()
