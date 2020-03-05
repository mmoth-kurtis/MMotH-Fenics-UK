---
page.title: Building a Mesh
---
One of the strengths of the finite element method is the ability to solve problems over complex geometries. Given cardiac MRI data, ventricle specific meshes can be generated via the following steps.

* Download MeVisLab version >= 3.1.1.  
* The following files from the MMotH-Fenics-UK repository are needed to generated meshes, and can be found in the [mesh_generation](https://github.com/mmoth-kurtis/MMotH-Fenics-UK/tree/master/mesh_generation) subdirectory:  
  * 1.createLV.py  
  * LV_Test.py  
  * LV_Segmentation_UKY.mlab  
* Open MeVisLab and load "LV_Segmentation_UKY.mlab". The screen should appear as follows:  
<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/mesh_image_01.png?raw=true" alt="titlepage" width="800"/>
* Using the DicomImport module in MeVisLab, convert raw MRI images to DICOM format.
* Follow instructions from FEniCS_LV_Segmentation.pdf. This should yield a '.stl' file. Do this for epi and endo volumes.
    * Some clarification for the Transform World Matrix: Change the entries in the last column so that the contours are in the center of the viewing volume. **Make sure to use the same transformation for both endo and epi volumes.**  
* Run the 1.createLV.py script on the stl files (creating json input for this now to be able to specify output directory and mesh name). This creates .vtk files.
* Run LV_Test.py on the New_Mesh files to yield the HDF5 file. This is what is needed for FEniCS. This step creates the mesh and assigns fiber angles, as well as the local coordinate system for each element.
* The mesh needs the top to be cropped out.



<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_01.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_02.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_03.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_04.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_05.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_06.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_07.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_08.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_09.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_10.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_11.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_12.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_13.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_14.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_15.jpg?raw=true" alt="titlepage" width="800"/>
