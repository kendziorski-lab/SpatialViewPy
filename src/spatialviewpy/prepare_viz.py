
def prepare10x_from_scanpy(adataObj,
                                   data_paths,
                                   export_path,
                                   samples_id_col = "sample_id",
                                   project_name = "spatial",
                                   cluster_col = "clusters",
                                   cluster_names = None,
                                   cluster_genes = None,
                                   clust_colors = None,
                                   layer = None,
                                   multisample_pattern = "_",
                                   expr_round = 3,
                                   spatial_sub_dir = "spatial",
                                   sample_info = None,
                                   download_repo = True,
                                   spatialview_repo = "https://github.com/kendziorski-lab/spatialview/archive/refs/tags/",
                                   spatialview_version = "spatialview-latest",
                                   port = 8878,
                                   launch_app = True,
                                   export_sparse = True,
                                   data_file_name_expressions = "expression_matrix.csv",
                                   data_file_name_expressions_sparse = "expression_matrix_sparse.txt",
                                   data_file_name_genes = "genes.csv",
                                   data_file_name_barcodes = "barcodes.csv",
                                   verbose = False,
                                   config_list = None):
    """
    Prepares data from Scanpy object and runs SpatialView.

    Args:
        adataObj (Scanpy): Scanpy object
        data_paths (str): directory path to raw data location
        export_path (str): directory path where the prepared data will be exported or SpatialView will be run
        samples_id_col (str, optional): Column name in adataObj containing sample names. Defaults to "sample_id".
        project_name (str, optional): Name of the project. Defaults to "spatial".
        cluster_col (str, optional): Column name in adataObj containing cluster annotations. Defaults to "clusters".
        cluster_names (str list, optional): Name of the clustes. If provided must match to the number of clusters. Defaults to None.
        cluster_genes (str list, optional): list of gene names with ',' separated. Defaults to None.
        clust_colors (str list, optional): list of hex colors corresponding to clusters. Defaults to None.
        layer (str, optional): Which later of Scanpy to be used. Defaults to None.
        multisample_pattern (str, optional): Character suffixed to the barcodes to make unique in case of multiple samples. Defaults to "_".
        expr_round (int, optional): The expression values are rounded upto. Defaults to 3.
        spatial_sub_dir (str, optional): subdirectory name where files related to spatial information are stored. Defaults to "spatial".
        sample_info (dataframe, optional): Metadata information about the samples. Rows corresponding to samples. Defaults to None.
        download_repo (bool, optional): If True, SpatialView will be downloaded from Github repository. Defaults to True.
        spatialview_repo (str, optional): SpatialView repo location. Defaults to "https://github.com/kendziorski-lab/spatialview/archive/refs/tags/".
        spatialview_version (str, optional): Which version to use. Defaults to "spatialview-latest".
        port (int, optional): Port number to be used by SpatialView. Defaults to 8878.
        launch_app (bool, optional): If True, then the SpatialView is launched. Defaults to True.
        export_sparse (bool, optional): Data compression for SpatialView. Defaults to True.
        data_file_name_expressions (str, optional): (For SpatialView) Name of the expression matrix. Defaults to "expression_matrix.csv".
        data_file_name_expressions_sparse (str, optional): (For SpatialView) Name of the expression matrix in sparse format. Defaults to "expression_matrix_sparse.txt".
        data_file_name_genes (str, optional): (For SpatialView) Name of the file containing genes name. Defaults to "genes.csv".
        data_file_name_barcodes (str, optional): (For SpatialView) Name of the file containing barcodes. Defaults to "barcodes.csv".
        verbose (bool, optional): Defaults to False.
        config_list (dictionary, optional): (For SpatialView) Additional configuration values to be passed to SpatialView. Defaults to None.

    Raises:
        FileNotFoundError: _description_
    """
    
    import numpy as np
    import pandas as pd
    from scipy.io import mmwrite
    from zipfile2 import ZipFile
    import seaborn as sns
    import requests
    import os
    import shutil
    import json
    import gzip
    from os import listdir, path
    
    orign_export_path = export_path
    samples_ids = adataObj.obs[samples_id_col].unique()
    sel_sample = None
    
    assert len(samples_ids) > 1 and len(samples_ids) == len(data_paths), "\
    Multiple samples found in the AnnData object. Length of data_paths is not same as number of samples.\n"
    
    # check if file path exists
    assert path.exists(export_path), f"{export_path} does not exist!"
    
    if download_repo:
        # TODO:: remove the url and provide the public url
        url = spatialview_repo + spatialview_version + '.zip'
        r = requests.get(url)
        
        if r.status_code == 200:
            with open(path.join(export_path,f'{project_name}.zip'), 'wb') as fh:
                fh.write(r.content)
        else:
            raise FileNotFoundError(f"Could not locate required SpatialView repository at {url}")

        #unzip the file
        tmp_zip = path.join(export_path, "spatialview_temp_")
        if path.exists(tmp_zip):
            shutil.rmtree(tmp_zip)
        os.mkdir(tmp_zip)
        ZipFile(path.join(export_path, f'{project_name}.zip')).extractall(path.join(export_path, "spatialview_temp_"))
        os.remove(path.join(export_path,f'{project_name}.zip'))
        spatialview_temp_dir = os.listdir(path.join(export_path, "spatialview_temp_"))
        #remove the existing project dir
        if path.exists(path.join(export_path, project_name)):
            shutil.rmtree(path.join(export_path, project_name))
        os.rename(path.join(export_path, "spatialview_temp_", spatialview_temp_dir[0]),
                 path.join(export_path, project_name))
        
        shutil.rmtree(path.join(export_path, "spatialview_temp_"))
        
        config_path = path.join(export_path, project_name, 'config', 'app_config.json')
        dataHTML_path = path.join(export_path, project_name, 'config', 'data_location.html')
        
        export_path = path.join(export_path, project_name, "data")
        
        spatalview_config = {}
        #updating the config file
        with open(config_path, 'r') as f:
            spatalview_config = json.load(f)
            if config_list :
                for k, v in config_list.items():
                    if k in spatalview_config:
                        spatalview_config[k] = v
                    
        spatalview_config['data_file_name_expressions'] = data_file_name_expressions
        spatalview_config['data_file_name_expressions_sparse'] = data_file_name_expressions_sparse
        spatalview_config['data_file_name_genes'] = data_file_name_genes
        spatalview_config['data_file_name_barcodes'] = data_file_name_barcodes
        spatalview_config['data_cluster_column'] = cluster_col
        
        with open(config_path, 'w') as json_file:
            json.dump(spatalview_config, json_file)
            
        #updating html data info tab
        dataInfo = pd.DataFrame({'Sample' : samples_ids})
        if sample_info is not None:
            dataInfo = pd.concat([dataInfo, dataInfo])
        else:
            sample_info = dataInfo
            
        dataInfo_html = dataInfo.to_html()
  
        # write html to file
        with open(dataHTML_path, "w") as f:
            f.write(dataInfo_html)
            
#     Processing each sample
    for i,sample_id in enumerate(samples_ids):
        sel_sample = sample_id
        data_path = data_paths[i]
        adataObj_temp = adataObj[adataObj.obs[samples_id_col] == sel_sample,]
        
        ## step 0
        #---------
        # create the export directory
        dir_name = sel_sample
        # remove space in names
        dir_name = dir_name.replace(" ", "")
        #first char should not not a number
        if dir_name[0].isdigit():
            dir_name = 'X' + dir_name
        
        #A name should be at least two characters long.
        if len(dir_name) == 1:
            dir_name = 'X' + dir_name
        
        if verbose and dir_name != sel_sample:
            print( f'{dir_name} is created for {sel.sel_sample}')
            
        export_path_sample = path.join(export_path, dir_name)
        os.mkdir(export_path_sample)
        
        # Step 1
        #--------
        # 1.a checking scalefactors_json.json, this file is a must have one else error
        sf_json_path = path.join(data_path, spatial_sub_dir, "scalefactors_json.json")
        sf_json_zip_path = path.join(data_path, spatial_sub_dir, "scalefactors_json.json.gz")
        
        assert path.exists(sf_json_path) or path.exists(sf_json_zip_path), f'scalefactors_json.json file not found at {sf_json_path}'
        
        if path.exists(sf_json_zip_path):
            with gzip.open(sf_json_zip_path,'rb') as f_in:
                with open(sf_json_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        shutil.copyfile(sf_json_path, path.join(export_path_sample, "scalefactors_json.json"))
        
        # 1.b checking tissue_positions_list.csv, this file is a must have one else error
        tp_csv_path = path.join(data_path, spatial_sub_dir, "tissue_positions_list.csv")
        tp_csv_zip_path = path.join(data_path, spatial_sub_dir, "tissue_positions_list.csv.gz")
        
        assert path.exists(tp_csv_path) or path.exists(tp_csv_zip_path), f'tissue_positions_list.csv file not found at {tp_csv_path}'
        
        if path.exists(tp_csv_zip_path):
            with gzip.open(tp_csv_zip_path,'rb') as f_in:
                with open(tp_csv_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        
        shutil.copyfile(tp_csv_path, path.join(export_path_sample, "tissue_positions_list.csv"))
        
        # Step 2
        #--------
        # Check either tissue_lowres_image.png or tissue_hires_image.png image present
        tlm_png = path.join(data_path, spatial_sub_dir, "tissue_lowres_image.png")
        thm_png = path.join(data_path, spatial_sub_dir, "tissue_hires_image.png")
        if not path.exists(thm_png) and not path.exists(tlm_png):
            raise f"Expected either {thm_png} or {tlm_png}"
        
        # Step 2.a
        # Checking if tissue_hires_image.png present or not, if not then trying to use tissue_lowres_image.png
        if not path.exists(thm_png) and path.exists(tlm_png):
            shutil.copyfileobj(tlm_png, thm_png)
            
            # updating the scalefactors_json.json file
            with open(sf_json_path) as sf_file:
                scalefactors = json.load(sf_file)
                scalefactors['tissue_hires_scalef'] = scalefactors['tissue_lowres_scalef']
                with open(path.join(export_path_sample, 'scalefactors_json.json'), "w") as outfile:
                    json.dump(scalefactors, outfile)
            
        shutil.copyfile(thm_png, path.join(export_path_sample, "tissue_hires_image.png"))
        
        # Step 3
        #--------
        # creating sample_info
        if sample_info is not None and isinstance(sample_info, pd.DataFrame):
            if sample_info.shape[0] >= i:
                sample_info_sample = pd.DataFrame(sample_info.loc[i,:])
                sample_info_sample.to_csv(path.join(export_path_sample, 'sample_info.csv'), 
                                        index = False, header = None)
            else:
                raise AssertionError("Number of rows in sample_info is not matching with number of samples")

        # Step 4
        #---------
        # creating metadata
        metadata = adataObj_temp.obs.copy()
        metadata['barcode'] = ["-".join(bc.split(multisample_pattern)[:-1]) for bc in metadata.index]
        
        if cluster_col in metadata.columns:
            metadata['cluster'] = metadata[cluster_col]
        else:
            metadata['cluster'] = 1
            
        spot_info = pd.read_csv(path.join(data_path, spatial_sub_dir,
                                          "tissue_positions_list.csv"),
                                header = None)        
        spot_info.columns = ["barcode", "in_tissue", "array_row", 
                              "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]  
        
        metadata = spot_info.merge(metadata, how='inner', on = "barcode")
        metadata.to_csv(path.join(export_path_sample, "metadata.csv"), index = False)
        
        # Step 5
        #--------
        #Extract the normalized counts
        
        #Note that we need (genes X cell) marix
        normalized_counts = adataObj_temp.X.T
        if layer:
            normalized_counts = adataObj_temp.layers[layer].T
        
        # creating expression_matrix data
        if export_sparse:
            normalized_counts = np.round(normalized_counts, expr_round)
            with open(path.join(export_path_sample,data_file_name_expressions_sparse), 'wb') as f:
                mmwrite(f,normalized_counts, precision = expr_round)
            
            # genes
            genes_df = pd.DataFrame(adataObj_temp.var_names)
            genes_df.to_csv(path.join(export_path_sample, data_file_name_genes), 
                            index = False, 
                            header= False)
            
            barcode_df = pd.DataFrame(["-".join(bc.split(multisample_pattern)[:-1]) for bc in adataObj_temp.obs.index])
            barcode_df.to_csv(path.join(export_path_sample, data_file_name_barcodes), 
                            index = False, 
                            header= False)
        else:
            normalized_counts = np.round(normalized_counts, expr_round)
            normalized_counts = pd.DataFrame(normalized_counts.toarray(), columns= adataObj_temp.var_names)
            normalized_counts['barcode'] = adataObj_temp.obs.index
            normalized_counts.to_csv(path.join(export_path_sample, data_file_name_expressions), index = False)
            
        # Step 6
        #----------
        # creating cluster info
        unique_clusters = 1
        if cluster_col in adataObj.obs.columns:
            unique_clusters = adataObj.obs[cluster_col].unique().sort_values()
        
        if not clust_colors:
            clust_colors = sns.color_palette("colorblind", len(unique_clusters))
            clust_colors = clust_colors.as_hex()
        
        if not cluster_names:
            cluster_names = unique_clusters

        if not cluster_genes:
            cluster_genes = ['undefined'] * len(unique_clusters)
            
        clusterInfo_df = pd.DataFrame({'cluster' : unique_clusters,
                                        'color' : clust_colors,
                                        'name' : cluster_names, 
                                        'genes' : cluster_genes})
        clusterInfo_df.to_csv(path.join(export_path_sample, 'cluster_info.csv'),
                              index = False)
            
    if download_repo and launch_app:
        start_httpserver(path.join(orign_export_path,project_name), port = port, verbose = False)
        

def start_httpserver(appPath, port = None, verbose = False, launch = True):
    from subprocess import Popen
    import webbrowser
    
    server = None
    http_port = 8000
    if not port:
        http_port = port
    try:
        server = Popen(['python', '-m', 'http.server', '--directory', appPath, str(http_port)])
        if launch:
            webbrowser.open_new_tab(f"http://localhost:{http_port}")
        if verbose:
            print(f'Process PID {server.pid}')
            print(f'Application is served at http://localhost:{http_port}')
    except Exception as e:
        print(e)
        
    return server