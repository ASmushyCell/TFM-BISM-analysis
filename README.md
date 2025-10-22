# READ ME FOR TRACTION FORCE MICROSCOPY (TFM) AND BAYESIAN INVERSION STRESS MICROSCOPY (BISM)

*This is a repository for the traction force microscopy and bayesian inversion stress microscopy pipeline used by the Ladoux-Mège group, Institut Jacques Monod, Paris.*

> [!NOTE]
> Written by Lucas Anger, building on the initial work of Philippe Marcq and Vincent Nier. Andreas Schönit and Fanny Wodrascka provided feedback and suggestions during the writing process. 


This code is an extension of the existing implementation of **Bayesian Inversion Stress Microscopy (BISM)**,  originally available at [pmarcq/BISM](https://github.com/pmarcq/BISM).
The main advantage of this pipeline is its user-friendly design. Assuming all dependencies are correctly installed, no prior coding experience is required to operate it. The main code performs the full analysis to compute both the traction force field and the stress field from beads displacement experiments. In addition, it provides basic, quick plots showing the temporal evolution of key parameters. Nonetheless, users with some experience in MATLAB or another programming language will find it easier to customize or troubleshoot the analysis when needed.

Relevant publications to read:

- [1] Dembo, M., Oliver, T., Ishihara, A., & Jacobson, K. (1996). Imaging the traction stresses exerted by locomoting cells with the elastic substratum method. Biophysical journal, 70(4), 2008-2022.
- [2] Trepat, X., Wasserman, M. R., Angelini, T. E., Millet, E., Weitz, D. A., Butler, J. P., & Fredberg, J. J. (2009). Physical forces during collective cell migration. Nature physics, 5(6), 426-430.
- [3] Tambe, D. T., Corey Hardin, C., Angelini, T. E., Rajendran, K., Park, C. Y., Serra-Picamal, X., ... & Trepat, X. (2011). Collective cell guidance by cooperative intercellular forces. Nature materials, 10(6), 469-475.
- [4] Nier, V., Jain, S., Lim, C. T., Ishihara, S., Ladoux, B., & Marcq, P. (2016). Inference of internal stress in a cell monolayer. Biophysical journal, 110(7), 1625-1635.
- [5] Saw, T. B., Doostmohammadi, A., Nier, V., Kocgozlu, L., Thampi, S., Toyama, Y., ... & Ladoux, B. (2017). Topological defects in epithelia govern cell death and extrusion. Nature, 544(7649), 212-216.
- [6] Teo, J. L., Lim, C. T., Yap, A. S., & Saw, T. B. (2020). A biologist’s guide to traction force microscopy using polydimethylsiloxane substrate for two-dimensional cell cultures. STAR protocols, 1(2), 100098.
- [7] Schoenit, A., Monfared, S., Anger, L., Rosse, C., Venkatesh, V., Balasubramaniam, L., ... & Ladoux, B. (2025). Force transmission is a master regulator of mechanical cell competition. Nature Materials, 1-11.

## CITATION

If you use this code, please cite the two following references:
- Nier, V., Jain, S., Lim, C. T., Ishihara, S., Ladoux, B., & Marcq, P. (2016). Inference of internal stress in a cell monolayer. Biophysical journal, 110(7), 1625-1635. [DOI:10.1016/j.bpj.2016.03.002](https://doi.org/10.1016/j.bpj.2016.03.002)
- Anger, L., Schoenit, A., Wodrascka, F., Rosse, C., Mège, R., Ladoux, B., Marcq, P. (2025). Tissue stress measurements with Bayesian Inversion Stress Microscopy, in revision, The European Physical Journal E.

## SOFTWARE REQUIREMENTS

The latest version of FIJI should be used, along with various plugins developed by external users — notably the Image Stabilizer and FTTC plugins.
For MATLAB, the analysis code was developed and tested using version R2024b, but it should be compatible with earlier versions, at least down to R2020. The PIVlab plugin must be installed in your MATLAB environment (preferably the latest version, although earlier versions should also work, at least from 2020 onwards).
Additionally, the pipeline makes use of several external MATLAB functions developed by other contributors. All required raw functions and FIJI plugins are already included in the repository’s functions folder. For reference, all original sources are summarized below.

### External sources :
  - **FIJI** 
    - **Image Stabilizer** plugin : K. Li, “The image stabilizer plugin for ImageJ,” [Image_Stabilizer](https://www.cs.cmu.edu/~kangli/code/Image_Stabilizer.html), February, 2008.
    - **FTTC** plugin : first developped for Tseng, Q. et al. Spatial Organization of the Extracellular Matrix Regulates Cell–cell Junction Positioning. PNAS (2012).doi:10.1073/pnas.1106377109, [LINK](https://sites.google.com/site/qingzongtseng/tfm)
  - **Matlab**
    - **PIVlab** plugin, open-source particle image velocimetry (PIV) software : first developped for Thielicke, W. and Stamhuis, E.J. (2014): PIVlab – Towards User-friendly, Affordable and Accurate Digital Particle Image Velocimetry in MATLAB. Journal of Open Research Software 2(1):e30, [LINK](https://www.PIVlab.de/) (you can also find it in the Matlab plugin library)
    - **200 colorbar**, user matlab package functions developped by Zhaoxu Liu / slandarer : [LINK](https://fr.mathworks.com/matlabcentral/fileexchange/120088-200-colormap)
    - **progressbar**, user matlab package functions developped by Steve Hoelzer : [LINK](https://fr.mathworks.com/matlabcentral/fileexchange/6922-progressbar)


## ABOUT THE EXEMPLE DATA 

A typical raw TFM dataset is provided in the `exemple_TFM` directory, located within the `data_test` folder. To generate this dataset, MDCK WT cells were seeded on a polydimethylsiloxane (PDMS) substrate with a stiffness of **15 kPa**, coated with fibronectin and 200 nm Cy3 fluorescent beads (FluoSpheres, Invitrogen), following the protocol fully described in [6]. \
The dataset includes the following image stacks:

- `brightfield` : a time-lapse stack of brightfield images focused on the cells.
- `WithCells` : a time-lapse stack of fluorescent bead images with cells on the substrate.
- `Ref` : a single reference image of the beads without cells (typically acquired at the end of the experiment, after removing the cells to allow the elastomer to relax).

Images were acquired every **15 minutes**, with a spatial resolution of **0.64 µm per pixel**.\
After completing the analysis, the contents of the `data_test` folder should match those of the `data_set_reference` folder, also located in `exemple_TFM`.

## STEPS OF THE PIPELINE

In order to use this pipeline, the raw data from a single position of a TFM experiment should be named according to the example data nomenclature:

- `brightfield` : is a time-lapse stack of brightfield images focused on the cells.
- `WithCells` : is a time-lapse stack of fluorescent bead images with cells on the substrate.
- `Ref` : is a single reference image of the beads without cells (typically acquired at the end of the experiment, after removing the cells to allow the elastomer to relax).

Provided your raw data are located in a specific folder, all required plugins are installed, and all functions are correctly linked to the matlab path, the pipeline can be run as follows:

- In **FIJI**, open the macro `pre_process_beforePIV.ijm` and run it. It will prompt you to select a folder; choose the folder containing the raw data you want to analyze. From here, the macro will run automatically. At the end of this step, a folder called `Input\forPIV` will be created, containing a sequence of images in which `WithCells` frames are interleaved with `Ref` frames for use in **PIVlab**.  
- In **MATLAB**, run the first section of the `tfm_main` script to reset the workspace. Ensure the current active folder is the one containing your raw data; if not, use the **Browse for folder** button.  
- From here, go to `Apps` and select **PIVlab_app**. The PIVlab interface will launch. Generally speaking, the interface is user-friendly, and users are encouraged to explore different parameters. We will detail the typical parameters used by the Ladoux-Mège team, but they can be changed depending on your experiments. In PIVlab, go to `File > New session`, then click the **Import images** button. In the window that opens, browse to the `Input\forPIV` folder of your working directory, select all images you want to analyze, and press **Add**. Once the full sequence appears in the _Images selected for import_ list, click **Import** to load them into PIVlab.
> [!IMPORTANT]
> During the import process, you can select an **Image sequencing style**. With this pipeline, the one you want to choose is the **_Pairwise_** sequencing style.
- At this point, your stack of images should occupy most of the screen. In the **Tools** box, you can navigate through each image. Go to `Analysis > PIV settings`: this is where you will set the parameters for the PIV analysis. If you want to explore the meaning of each parameter, we encourage you to consult the PIVlab documentation. Typically, in our TFM experiments, we only adjust the number of passes and the size of the corresponding interrogation areas. For the example dataset, we recommend using a first window of **64 pixels** with an overlap of **32 pixels**, and a second window of **32 pixels** with an overlap of **16 pixels** (note that the optimal settings may vary depending on your microscope camera resolution and magnification).  This will produce an output field with one value for every 16x16 pixel square in the original image. You can preview the analysis on the current active frame to fine-tune parameters by clicking the **Analyse current frame** button. Once satisfied, go to `Analysis > ANALYSE!` and run the **Analyse all frames** button.
- Next, we filter and improve the coherence of the computed displacement field, as some values may be aberrant. Go to `Validation > Velocity based validation`. A panel will open with various parameters for smoothing. In our lab, we typically leave the default statistical smoothing parameters unchanged. In addition to statistical filtering, there is a manual filter based on a brute-force approach: click **Refine velocity limits** and drag the cursor to create a rectangular window around the cloud of points representing all displacements in the x and y directions. The goal is to include the global cloud while excluding clear outliers. Once done, click **Apply to all frames**.
- The PIV analysis is now complete. You have obtained the displacement field of the beads for your current folder. Save the result by going to `File > Export > MAT file`, click **Export all frames**, select your current working folder, and save the file under the name **PIVlab.mat** (this should be the default name). You can then close PIVlab and return to `tfm_main.mat`.
- From this point, run each section of the script sequentially. After completing the BISM section, the remaining sections are primarily dedicated to data visualization, allowing you to plot various parameters over time and space. These sections can be skipped depending on your analysis goals and needs.

## SUMMARY OF OUTPUTS

> [!NOTE]
> - Each cell in the following variables corresponds to one frame of your initial image stack, where **N** is the total number of frames.  
> - Field dimensions are defined by the **PIV grid** generated in the *PIVlab* plugin steps, both by the window size and the overlap.  
> - **stress_yx** is not computed, as by construction, the stress tensor is **symmetric**.
> - **stress_xy** is also referred to as **shear stress**. Since it is not invariant under a change of basis, if you wish to quantify the overall "amount" of shear stress, you should instead use an invariant quantity such as **stress_aniso**, sometimes also called the **maximum shear stress**.
> - Boundary values are set to `NaN` to remove unreliable estimates at the image edges.
> - This pipeline assumes that the tissue behaves as a 2D sheet. As a result, some of the computed parameters may appear to have inconsistent dimensions from a physicist’s perspective, notably including an extra unit length [L] factor. This is because they are calculated under the assumption of a constant and known tissue height. For a detailed explanation, please refer to [4]. If you wish to express these variables in their proper physical units, you can divide them by the tissue height, provided it remains reasonably uniform in magnitude and its value is known.

- **`traction_force.mat`**: This file contains the **traction force field** computed from FTTC (Fourier Transform Traction Cytometry) analysis.  
It includes the following variables, vector components or derived quantities for each time frame of the experiment:

| Variable | Description | Type |
|-----------|--------------|------|
| **Tx** | **x-component** of the traction force field, in pascal (Pa) | Cell array of 2D matrices, of size {N,1} |
| **Ty** | **y-component** of the traction force field, in pascal (Pa) | Cell array of 2D matrices, of size {N,1} |
| **T** | Traction force **magnitude** field, in pascal (Pa). | Cell array of 2D matrices, of size {N,1} |
| **x_TFM** | **x-coordinates** in pixels in the original image corresponding to each traction value. | Cell array of 2D matrices, of size {N,1} |
| **y_TFM** | **y-coordinates** in pixels in the original image corresponding to each traction value. | Cell array of 2D matrices, of size {N,1} |

- **`stress_mean.mat`**: This file contains the **stress tensor field** computed from the traction force field using BISM.  
It includes the following variables, corresponding to the reconstructed stress tensor components and derived scalar quantities for each time frame of the experiment:

| Variable | Description | Type |
|-----------|--------------|------|
| **stress_xx** | **xx-component** of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_yy** | **yy-component** of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_xy** | **xy-component** of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_M** | **Maximum principal stress**, the larger eigenvalue of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_m** | **Minimum principal stress**, the smaller eigenvalue of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_iso** | **Isotropic stress**, first invariant of the stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **stress_aniso** | **Anisotropic stress**, second invariant of the deviatoric stress tensor, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **sigma_VM** | **Von Mises equivalent stress**, in pascal micron (Pa.µm). | Cell array of 2D matrices, of size {N,1} |
| **angle_stress_x**, **angle_stress_y** | unit x-component and y-component of the **principal stress orientation** director field. | Cell arrays of 2D matrices, of size {N,1} |

- **`parameter_stats.mat`**: This file contains time-dependent statistical characterization of the parameters computed using TFM and BISM analyses.
It includes the following variables:
  
| Variable | Description | Type |
|-----------|--------------|------|
| **time** | Time points for each frame, in hour (h). | Array of double of size 1 x N |
| **m_T** | Mean **traction force magnitude** across the field for each time frame, in pascal (Pa). | Array of double of size N x 1 |
| **f_T** | Standard deviation of **traction force magnitude** across the field for each time frame, in pascal (Pa). | Array of double of size N x 1 |
| **m_iso** | Mean **isotropic stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |
| **f_iso** | Standard deviation of **isotropic stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |
| **m_aniso** | Mean **anisotropic stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |
| **f_aniso** | Standard deviation of **anisotropic stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |
| **m_vm** | Mean **von Mises stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |
| **f_vm** | Standard deviation of **von Mises stress** across the field for each time frame, in pascal micron (Pa.µm). | Array of double of size N x 1 |

## IMPORTANT NOTES, TROUBLESHOOTING AND COMMON ISSUES

- In some cases, the bead images may appear slightly shaky, with small shifts in X and Y, as microscopes are not always perfectly stable. In such situations, a stronger alignment plugin for FIJI can help correct these drifts. You can find it here: [Linear Stack Alignment with SIFT](https://imagej.net/plugins/linear-stack-alignment-with-sift).  It is recommended to use this plugin before running the rest of the pipeline.


> [!TIP]
> Organize different datasets by placing each one in a separate directory.  
> Typically, for a given experiment, the structure should initally look like this:
>
> ```
> name_of_your_experiment/
> ├── pos1/
> │   ├── brightfield.tif
> │   ├── Ref.tif
> │   └── WithCells.tif
> ├── pos2/
> │   ├── brightfield.tif
> │   ├── Ref.tif
> │   └── WithCells.tif
> └── ...
> ```
> Do **not** rename the output files generated by the different scripts unless you are confident working with MATLAB. Many scripts assume specific filenames for proper execution, and keeping filenames consistent also allows you to run batch analyses efficiently across multiple datasets!

- Most scripts in this repository contain additional commented lines, often including notes that describe the purpose of each section and variable. If you encounter errors or unexpected results, review the commented sections : they may help you diagnose the issue. As a rule of thumb, carefully **read the error messages**: they often provide useful hints about the source of the problem.
- A common source of errors is missing dependencies, such as:
   - Functions not added to the current MATLAB path.  
   - Fiji plugins installed in the wrong plugin directory.  
- Another frequent issue occurs during the initial analysis step: the **ImageJ “Invert LUT”** option may have been applied to the `Ref` frame but not to `WithCells` frames, or vice versa.  This mismatch can lead to incorrect displacement maps.
- Even if everything works as expected, you should always ask yourself: do the values that were just computed make sense? Do they have a physical meaning that seems reasonable?

  




















