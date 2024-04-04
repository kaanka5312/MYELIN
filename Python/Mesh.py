import nibabel as nib
import numpy as np

# Load the parcellation labels
parcellation_path = 'E:/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii'
parcellation_data = nib.load(parcellation_path)

data, brain_models = parcellation_data.get_fdata(), parcellation_data.header.get_index_map(1).brain_models

for bm in brain_models:
    if bm.model_type == 'CIFTI_MODEL_TYPE_SURFACE':
        # This is a surface model, could be cortex left, cortex right
        print(f"Model Type: {bm.model_type}, Brain Structure: {bm.brain_structure}")
        # Now you can access the vertex indices for this brain structure
        vertex_indices = bm.vertex_indices
        # And extract the data for these vertices
        surface_data = data[vertex_indices]



# Load the surface mesh
surface_path = 'your_surface_file.gii'
surface_mesh = nib.load(surface_path)
vertices = surface_mesh.darrays[0].data
faces = surface_mesh.darrays[1].data

# Specify the parcel ID you're interested in extracting
parcel_id = 1  # Example parcel ID

# Find vertices that belong to the specified parcel
parcel_vertices_indices = np.where(parcellation_data == parcel_id)[0]

# Extract the faces that only contain vertices within the specified parcel
parcel_faces = np.array([face for face in faces if np.all(np.isin(face, parcel_vertices_indices))])

# Adjust the indices of the vertices to start from 0 for the new sub-mesh
unique_vertices_indices = np.unique(parcel_faces)
new_indices_map = {old_idx: new_idx for new_idx, old_idx in enumerate(unique_vertices_indices)}
new_faces = np.vectorize(new_indices_map.get)(parcel_faces)

# Extract the corresponding vertices
new_vertices = vertices[unique_vertices_indices]

# At this point, `new_vertices` and `new_faces` define the sub-mesh for the specified parcel.
# You can then save these as a new GIFTI file or use them directly in your analysis.
