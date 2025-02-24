import numpy as np

# A sandbox script to define an FFD box and rotate it about an arbitrary center of rotation and around a unit vector
# using Rodrigues rotation formula

# Basic steps involved:  1) Translate the vertices so that center of rotation becomes the origin
#                        2) Perform the rotation around a unit vector (For example Y axis ) 
#                        3) Translate the vertices back to their original position relative to the specified center


def ffd_box(vertices):

    # Define the faces of the parallelepiped (using vertex indices)
    faces = [
        [0, 1, 2, 3],  # Lower surface
        [4, 5, 6, 7],  # Upper surface
        [0, 1, 5, 4],  # Side 1
        [1, 2, 6, 5],  # Side 2
        [2, 3, 7, 6],  # Side 3
        [3, 0, 4, 7]   # Side 4
    ]
    return faces
#end

def write_ffd_box(faces, vertices, filename):
    # Write the data to a .dat file in Tecplot ASCII format
    with open(filename, 'w') as f:
        # Header
        f.write("TITLE = \"Parallelepiped Geometry\"\n")
        f.write("VARIABLES = \"X\", \"Y\", \"Z\"\n")
        f.write(f"ZONE N={len(vertices)}, E={len(faces)}, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n")
    
        # Write vertices
        for vertex in vertices:
            f.write(f"{vertex[0]} {vertex[1]} {vertex[2]}\n")
    
        # Write faces (connectivity, using 1-based indexing for Tecplot)
        for face in faces:
            f.write(f"{face[0]+1} {face[1]+1} {face[2]+1} {face[3]+1}\n")

    print(f"Parallelepiped data has been written to {filename}.")
#end

def rotate_vertices(vertices, axis, theta_degrees, center):
    """Rotate vertices around an arbitrary axis vector by an angle theta (in radians)."""
    # Ensure the axis is a unit vector
    axis = axis / np.linalg.norm(axis)

    print("Rotation axis: ", axis)
    
    # Calculate rotation matrix using Rodrigues' rotation formula
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])

    # Translate vertices so that center of rotation becomes the origin.
    translated_vertices = vertices - center

    theta =np.radians(theta_degrees)
    
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * np.dot(K, K)
    
    # Rotate all vertices
    rotated_translated_vertices = np.dot(translated_vertices, R.T)

    # Translate vertices back to the original position
    rotated_vertices = rotated_translated_vertices + center
    return rotated_vertices


def compute_face_normal(v0, v1, v2, v3):
    """
    Compute the unit normal vector of a face given four vertices.

    Parameters:
    v0, v1, v2, v3 : array-like
        The coordinates of the four vertices defining the face. Each should be a 3-element array-like structure.

    Returns:
    unit_normal : ndarray
        The unit normal vector of the face.
    """
    # Compute two vectors from the three vertices
    vector1 = np.array(v1) - np.array(v0)
    vector2 = np.array(v2) - np.array(v0)
    
    # Calculate the normal using the cross product
    normal = np.cross(vector1, vector2)
    
    # Normalize the normal vector to get the unit normal
    unit_normal = normal / np.linalg.norm(normal)
    return unit_normal

# Define the vertices of the parallelepiped in counter-clockwise order starting from the lower surface

vertices = np.array([
    [0.0, 0.3, 0.0],  # Lower root LE
    [1.0, 0.3, 0.0],  # Lower root TE
    [1.0, 4.0, 0.5],  # Lower tip TE
    [0.0, 4.0, 0.5],  # Lower tip LE
    [0.0, 0.0, 1.0],  # Upper root LE
    [1.0, 0.0, 1.0],  # Upper root TE
    [1.0, 4.0, 1.5],  # Upper tip TE
    [0.0, 4.0, 1.5]   # Upper tip LE
], dtype=np.float32)

# Generate the FFD box
faces = ffd_box(vertices)

# Write the original FFD box
write_ffd_box(faces, vertices, 'ffd_box_0.dat')

# Compute the face normals
v0 = [0.0, 0.3, 0.0]
v1 = [1.0, 0.3, 0.0]
v2 = [0.0, 0.0, 1.0]
v3 = [1.0, 0.0, 1.0]


# Compute the unit normal vector for the face
unit_normal = compute_face_normal(v0, v1, v2, v3)


theta = 5

#rotation_vector = np.array([0.0, 1.0, 0.0])  # Rotate around the Y-axis
# The Rotate the FFD about about the surface normal of the root face
rotation_vector = unit_normal
center_of_rotation = np.array([0.5, 0.0, 0.5])  # Center of rotation



# Rotate the vertices
rotated_vertices = rotate_vertices(vertices, rotation_vector, theta, center_of_rotation)

# Write the rotated FFD box
write_ffd_box(faces, rotated_vertices, 'ffd_box_rotated.dat')
