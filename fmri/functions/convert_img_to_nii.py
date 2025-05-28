import nibabel as nb

def convert_img_to_nii(path_to_file):
    """
    Converts a .img file to a .nii file
    """
    img = nb.load(path_to_file)
    nb.save(img, path_to_file.replace('.img', '.nii'))
