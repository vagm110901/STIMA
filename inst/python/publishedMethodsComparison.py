# -*- coding: utf-8 -*-

def pythonEnviroment(mode):
    """
    Import the python environment information.

    Parameters
    ----------
    mode:  str. A value in ["STalign", "PASTE2"]
    """

    if mode not in ["STalign", "PASTE2"]:
        raise ValueError("mode value must be 'STalign' or 'PASTE2'.")

    if mode == "STalign":
        import sys
        import numpy as np
        import STalign
        import matplotlib.pyplot as plt 
        import torch
    elif mode == "PASTE2":
        import os
        import pandas as pd
        import numpy as np
        import scanpy as sc
        import seaborn as sns
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from scipy.sparse import csr_matrix
        from skimage import io
        from skimage.draw import disk
        import anndata as ad
        
        
        
    

def normalize_images(image):
    """
    Normalize the image to the range [0, 1] using STalign function to use it in R.
    
    Parameters
    ----------
    image:  torch tensor
            Image to be normalized.
    
    Returns
    -------
    im_norm: torch tensor
             Normalized image.
    """
    im_norm = STalign.normalize(image)
    return im_norm


def STalign_transformation(points_im_prob_rc, points_im_ref_rc, 
                           im_prob_norm, im_ref_norm
                           y_im_prob, x_im_prob):
    """
    Perform the STalign transformation on the images and points.
    
    Parameters
    ----------
    points_im_prob_rc: torch tensor
                       N x 2 set of corresponding points for matching in source image.
    points_im_ref_rc:  torch tensor
                       N x 2 set of corresponding points for matching in target image.
    im_prob_norm: torch tensor
                  Normalized source image.
    im_ref_norm:  torch tensor
                  Normalized target image.
    y_im_prob: torch tensor
               Y coordinates of the source image.
    x_im_prob: torch tensor
               X coordinates of the source image.

    Returns
    -------
    align_image_im_prob:  torch tensor
                          Aligned source image.
    align_points_im_prob: torch tensor
                          Aligned points in the source image. 
    """

    # compute initial affine transformation from points
    L,T = STalign.L_T_from_points(points_im_prob_rc, points_im_ref_rc)

    # transpose matrices into right dimensions, set up extent 
    im_ref_trans = im_ref_norm.transpose(2,0,1)
    Y_im_ref = np.array(range(im_ref_trans.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    X_im_ref = np.array(range(im_ref_trans.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    extent_im_ref = STalign.extent_from_x((Y_im_ref,X_im_ref))

    # transpose matrices into right dimensions, set up extent
    im_prob_trans = im_prob_norm.transpose(2,0,1)
    Y_im_prob = np.array(range(im_prob_trans.shape[1]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    X_im_prob = np.array(range(im_prob_trans.shape[2]))*1. # needs to be longs not doubles for STalign.transform later so multiply by 1.
    extent_im_prob = STalign.extent_from_x((Y_im_prob,X_im_prob))

    torch.set_default_device('cpu')
    device = 'cpu'

    # keep all other parameters default
    params = {'L':L,'T':T,
        'niter':400,
        'pointsI':points_im_prob_rc,
        'pointsJ':points_im_ref_rc,
        'device':device,
        'sigmaM':0.15,
        'sigmaB':0.2,
        'sigmaA':0.3,
        'a':300,
        'epV':10,
        'muB': torch.tensor([255,255,255]) # white is background
    }

    # run LDDMM
    out = STalign.LDDMM([Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref],im_ref_trans,**params)

    # get necessary output variables
    A = out['A']
    v = out['v']
    xv = out['xv']

    # modify the source
    align_points_im_prob = STalign.transform_points_source_to_target(xv,v,A,np.stack([y_im_prob, x_im_prob], -1)) # transform the points
    align_points_im_prob = align_points_im_prob.numpy() # convert from tensor to numpy to access in R later

    im_prob_trans = im_prob_trans.astype(np.float64)
    align_image_im_prob = STalign.transform_image_source_to_target(xv,v,A,[Y_im_prob,X_im_prob],im_prob_trans,[Y_im_ref,X_im_ref])
    align_image_im_prob = align_image_im_prob.permute(1,2,0).numpy()

    return {
        "align_image_im_prob": align_image_im_prob,
        "align_points_im_prob": align_points_im_prob}


def h5ad_create():
    """
    Create an AnnData (h5ad) file using the r.saveSeurat_forAnnData() function returns.
    It is run in Python, and it finds the necessary files in "./results/AnnData/".
    The return is a file saved in the same directory.
    
    Notes:
        - The function reads from the directory "./results/AnnData/".
        - It saves the h5ad AnnData files to "./results/AnnData/".
        - It relies on the naming convention "expression_matrix_slice{i}.csv", "cell_metadata_slice{i}.csv",
          "gene_metadata{i}.csv", "image_coordinates_slice{i}.csv", "scale_factors_slice{i}.csv",
          "spatial_image_slice{i}.png" and "PASTE2_merge_{i}.h5ad".
    """
    pythonEnviroment("PASTE2")

    carpetaData = "./results/AnnData/"

    slice_indices = sorted(set(f.split("slice")[1].split(".")[0] for f in os.listdir(saveDir) if "slice" in f))
    for i in slice_indices:
        i = str(i)
        # Load expression matrixes and metadata for each slice
        expression = pd.read_csv(f'{carpetaData}expression_matrix_slice{i}.csv', index_col=0)
        cell_metadata = pd.read_csv(f'{carpetaData}cell_metadata_slice{i}.csv', index_col=0)
        gene_metadata = pd.read_csv(f'{carpetaData}gene_metadata{i}.csv', index_col=0)  
        
        # Cargar las coordenadas espaciales (row, col, imagerow, imagecol)
        image_coords = pd.read_csv(f'{carpetaData}image_coordinates_slice{i}.csv', index_col=0)

        image_coords[["col", "row"]] = image_coords[["row", "col"]]
        image_coords[["imagerow", "imagecol"]] = image_coords[["imagecol", "imagerow"]]

        spatial_coords = image_coords[['imagerow', 'imagecol']]  # Ajustar según los nombres de las columnas
        scale_factors = pd.read_csv(f'{carpetaData}scale_factors_slice{i}.csv', index_col=0)
        
        image_coords_selected = image_coords.iloc[:,:]
        cell_metadata = cell_metadata.join(image_coords_selected, how="left")
        
        # Cargar la imagen espacial para el slice
        spatial_image = io.imread(f'{carpetaData}spatial_image_slice{i}.png')

        # Crear objeto AnnData para el slice
        adata = sc.AnnData(X=expression.values, obs=cell_metadata, var=gene_metadata)

        # Incluir las coordenadas espaciales
        adata.obsm['spatial'] = spatial_coords.values
        #adata.obsm['rgb'] = spatial_image

        # Asociar la imagen al objeto AnnData
        image_name = str(cell_metadata.name.iloc[0])
        adata.uns['spatial'] = {image_name: {}}
        adata.uns['spatial'][image_name]['images'] = {}
        adata.uns['spatial'][image_name]['images'] = {'lowres': spatial_image}
        adata.uns['spatial'][image_name]['scalefactors'] = {
            'tissue_hires_scalef': float(scale_factors.tissue_hires_scalef.iloc[0]),
            'tissue_lowres_scalef': float(scale_factors.tissue_lowres_scalef.iloc[0]),
            'spot_diameter_fullres': float(scale_factors.spot_diameter_fullres.iloc[0]),
            'fiducial_diameter_fullres': float(scale_factors.fiducial_diameter_fullres.iloc[0])}
        
        adata.var_names = adata.var['x']
        adata.X = csr_matrix(adata.X)

        adata.write(f'{carpetaData}PASTE2_merge_{i}.h5ad')


def RGBvalues():
    """
    Extracts and assigns RGB values to spatial transcriptomics data.

    This function reads `.h5ad` AnnData files and corresponding low-resolution spatial images from the specified
    directory, calculates the average RGB values for each spot (cell) based on their spatial coordinates, and 
    saves the updated AnnData objects with the extracted RGB data.

    Notes:
        - The function reads from the directory "./results/AnnData/".
        - It saves the updated h5ad AnnData files to "./results/PASTE2/".
        - It relies on the naming convention "PASTE2_merge_{i}.h5ad" and "spatial_image_slice{i}.png".
    """
    pythonEnviroment("PASTE2")
    import cv2

    carpetaData = "./results/AnnData/"
    saveDir = "./results/PASTE2/"

    adatas = []
    ims = []
    slice_indices = sorted(set(f.split("slice")[1].split(".")[0] for f in os.listdir(carpetaData) if "slice" in f))
    for i in slice_indices:
        adata_name = f"{carpetaData}PASTE2_merge_{i}.h5ad"
        adata = sc.read_h5ad(adata_name)
        adatas.append(adata)
        im_name = f"{carpetaData}spatial_image_slice{i}.png"
        im = cv2.imread(im_name)
        ims.append(im)

    for i in range(0,len(adatas)):
        adata = adatas[i]
        im = ims[i]
        # I am using the lowres image, because is the one that was available 
        image_name = str(adata.obs.name.iloc[0])
        spot_diam_fullres = float(adata.uns['spatial'][image_name]['scalefactors']['spot_diameter_fullres'])
        scale_lowres = float(adata.uns['spatial'][image_name]['scalefactors']['tissue_lowres_scalef'])
        spot_diam_lowres = round(spot_diam_fullres * scale_lowres)
        spot_rad_lowres = round(spot_diam_lowres / 2)

        rgb = []
        for j in range(adata.n_obs):
            x, y = round(adata.obsm['spatial'][j][0] * scale_lowres), round(adata.obsm['spatial'][j][1] * scale_lowres)
            spot_mask = np.zeros((im.shape[0], im.shape[1]), np.uint8)
            cv2.circle(spot_mask, (x, y), radius=int(spot_rad_lowres), color=(255, 255, 255), thickness=-1)
            rgb.append(tuple(int(round(c)) for c in cv2.mean(im, spot_mask)[::-1][1:]))
            
        adata.obsm['rgb'] = np.array(rgb)

        adata.write(f'{saveDir}PASTE2_merge_{i+1}.h5ad')


def partial_pairwise_align_histology(sliceA, sliceB, alpha=0.1, s=None, armijo=False, dissimilarity='glmpca', use_rep=None, G_init=None, a_distribution=None,
                   b_distribution=None, norm=True, return_obj=False, verbose=False, **kwargs):
    """
    Optimal partial alignment of two slices using both gene expression and histological image information.

    sliceA, sliceB must be AnnData objects that contain .obsm['rgb'], which stores the RGB value of each spot in the histology image.
    """
    m = s
    print("PASTE2 starts...")

    # subset for common genes
    common_genes = intersect(sliceA.var.index, sliceB.var.index)
    sliceA = sliceA[:, common_genes]
    sliceB = sliceB[:, common_genes]
    # print('Filtered all slices for common genes. There are ' + str(len(common_genes)) + ' common genes.')

    # Calculate spatial distances
    D_A = distance.cdist(sliceA.obsm['spatial'], sliceA.obsm['spatial'])
    D_B = distance.cdist(sliceB.obsm['spatial'], sliceB.obsm['spatial'])

    # Calculate expression dissimilarity
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, use_rep)), to_dense_array(extract_data_matrix(sliceB, use_rep))
    if dissimilarity.lower() == 'euclidean' or dissimilarity.lower() == 'euc':
        M_exp = distance.cdist(A_X, B_X)
    elif dissimilarity.lower() == 'kl':
        s_A = A_X + 0.01
        s_B = B_X + 0.01
        M_exp = kl_divergence(s_A, s_B)
    elif dissimilarity.lower() == 'glmpca':
        M_exp = glmpca_distance(A_X, B_X, latent_dim=50, filter=True, verbose=verbose)
    else:
        print("ERROR")
        exit(1)

    # Calculate RGB dissimilarity
    # sliceA_rgb = (sliceA.obsm['rgb'] - np.mean(sliceA.obsm['rgb'], axis=0)) / np.std(sliceA.obsm['rgb'], axis=0)
    # sliceB_rgb = (sliceB.obsm['rgb'] - np.mean(sliceB.obsm['rgb'], axis=0)) / np.std(sliceB.obsm['rgb'], axis=0)
    M_rgb = distance.cdist(sliceA.obsm['rgb'], sliceB.obsm['rgb'])
    # M_rgb = distance.cdist(sliceA_rgb, sliceB_rgb)

    # Scale M_exp and M_rgb, obtain M by taking half from each
    M_rgb /= M_rgb[M_rgb > 0].max()
    M_rgb *= M_exp.max()
    # M_exp /= M_exp[M_exp > 0].max()
    # M_rgb /= M_rgb[M_rgb > 0].max()
    
    M = 0 * M_exp + 1 * M_rgb               ### line 345
    # M = 0.5 * M_exp + 0.5 * M_rgb 

    # init distributions
    if a_distribution is None:
        a = np.ones((sliceA.shape[0],)) / sliceA.shape[0]
    else:
        a = a_distribution

    if b_distribution is None:
        b = np.ones((sliceB.shape[0],)) / sliceB.shape[0]
    else:
        b = b_distribution

    if norm:
        D_A /= D_A[D_A > 0].min().min()
        D_B /= D_B[D_B > 0].min().min()

        """
        Code for normalizing distance matrix
        """
        D_A /= D_A[D_A>0].max()
        D_A *= M.max()
        D_B /= D_B[D_B>0].max()
        D_B *= M.max()
        """
        Code for normalizing distance matrix ends
        """

    # Run OT
    pi, log = PASTE2.partial_fused_gromov_wasserstein(M, D_A, D_B, a, b, alpha=alpha, m=m, G0=G_init, loss_fun='square_loss', armijo=armijo, log=True, verbose=verbose)

    if return_obj:
        return pi, log['partial_fgw_cost']
    return pi


def PASTE2_align():
    """
    Performs the alignemnt by PASTE2.

    Notes:
        - The function reads from the directory "./results/PASTE2/".
        - It saves the updated h5ad AnnData PASTE2-aligned files to "./results/PASTE2/".
        - It relies on the naming convention "PASTE2_merge_{i}.h5ad" and "PASTE2_merge_{i}_align1{i}_h.h5ad".
    """
    pythonEnviroment("PASTE2")
    from paste2 import PASTE2, projection, model_selection
    from scipy.spatial import distance
    from paste2.helper import kl_divergence, intersect, to_dense_array, extract_data_matrix, generalized_kl_divergence, high_umi_gene_distance, pca_distance, glmpca_distance

    carpetaData = "./results/PASTE2/"

    slice_indices = sorted(set(f.split("align1")[1].split("_h.")[0] for f in os.listdir(carpetaData) if "slice" in f))

    patient_1_name = f"{carpetaData}PASTE2_merge_1.h5ad"
    patient_1 = sc.read_h5ad(patient_1_name)

    rowmin_0 = min(patient_1.obs.imagerow)
    colmin_0 = min(patient_1.obs.imagecol)
    rowmax_0 = max(patient_1.obs.imagerow)
    colmax_0 = max(patient_1.obs.imagecol)

    for i in slice_indices:
        if i != 1:
            patient_prob_name = f"{carpetaData}PASTE2_merge_{i}.h5ad"
            patient_prob = sc.read_h5ad(patient_prob_name)

            rowmax = max(patient_prob.obs.imagerow)
            colmax = max(patient_prob.obs.imagecol)
            patient_prob.obsm['spatial'][:,0] = patient_prob.obsm['spatial'][:,0] / rowmax * rowmax_0
            patient_prob.obsm['spatial'][:,1] = patient_prob.obsm['spatial'][:,1] / colmax * colmax_0

            s_1prob = 0.9

            # 100% Histology
            pi_h = partial_pairwise_align_histology(patient_1, patient_prob, s = s_1prob)
            new_h = projection.partial_stack_slices_pairwise([patient_1, patient_prob], [pi_h])
            patient_1_align = new_h[0]
            patient_prob_align = new_h[1]

            rowmin = np.min(patient_1_align.obsm['spatial'][:,0])
            colmin = np.min(patient_1_align.obsm['spatial'][:,1])
            patient_1_align.obsm['spatial'][:,0] = patient_1_align.obsm['spatial'][:,0] - rowmin + rowmin_0
            patient_1_align.obsm['spatial'][:,1] = patient_1_align.obsm['spatial'][:,1] - colmin + colmin_0
            patient_prob_align.obsm['spatial'][:,0] = patient_prob_align.obsm['spatial'][:,0] - rowmin + rowmin_0
            patient_prob_align.obsm['spatial'][:,1] = patient_prob_align.obsm['spatial'][:,1] - colmin + colmin_0

            patient_1_align.write(f"{saveData}PASTE2_merge_1_align1{N}_h.h5ad")
            patient_prob_align.write(f"{saveData}PASTE2_merge_{N}_align1{N}_h.h5ad")
    

def saveAnnData_forSeurat():
    """
    Save csv files with the coordinates using AnnData files after alignment. 

    Notes:
        - The function reads from the directory "./results/PASTE2/".
        - It saves the csv files to "./results/PASTE2/".
        - It relies on the naming convention "PASTE2_merge_{i}_align1{i}_h.h5ad" and "PASTE2_{i}_align1{i}_h_coord.csv".
    """
    pythonEnviroment("PASTE2")

    carpetaData = "./results/PASTE2/"

    slice_indices = sorted(set(f.split("align1")[1].split("_h.")[0] for f in os.listdir(carpetaData) if "slice" in f))
    for i in slice_indices:
        patient_1_name = (f"{carpetaData}PASTE2_merge_1_align1{i}_h.h5ad")
        patient_1 = sc.read_h5ad(patient_1_name)

        patient_prob_name = (f"{carpetaData}PASTE2_merge_{i}_align1{i}_h.h5ad")
        patient_prob = sc.read_h5ad(patient_prob_name)

        # 2 columns: imagerow, imagecol
        np.savetxt(f"{carpetaData}PASTE2_1_align1{i}_h_coord.csv",
                patient_1.obsm['spatial'], delimiter=",", fmt="%.2f", header="imagerow,imagecol", comments="")
        np.savetxt(f"{carpetaData}PASTE2_{i}_align1{i}_h_coord.csv",
                patient_prob.obsm['spatial'], delimiter=",", fmt="%.2f", header="imagerow,imagecol", comments="")
