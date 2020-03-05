---
page.title: Building a Mesh
---
One of the strengths of the finite element method is the ability to solve problems over complex geometries. Given cardiac MRI (CMR) data, ventricle specific meshes can be generated.

First, download [MeVisLab](https://www.mevislab.de/download/) version >= 3.1.1. Any previous version will be lacking required modules. Additionally, the following files are required from the MMotH-Fenics-UK repository and can be found in the [mesh_generation](https://github.com/mmoth-kurtis/MMotH-Fenics-UK/tree/master/mesh_generation) subdirectory:  
  * 1.createLV.py  
  * LV_Test.py  
  * LV_Segmentation_UKY.mlab  

Open MeVisLab and load "LV_Segmentation_UKY.mlab". The following screen should appear:

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/mesh_image_02.png?raw=true" alt="titlepage" width="800"/>  

Go ahead and exit out of the "Panel WEMSave" box if it appears. MeVisLab saves the state the last time this file was opened, and the WEMSave module is the last to be called before exiting.  

First, we need to convert the raw CMR data to the DICOM format:  

* Convert the raw CMR data to the DICOM format:  
  * Double click the DicomImport module to open it. Navigate to the directory containing the raw CMR images.
  * Specify the folder to save the converted images in the "Target Path" dialogue box.
* Follow instructions from FEniCS_LV_Segmentation.pdf. This should yield a '.stl' file. Do this for epi and endo volumes.

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_03.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_04.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_05.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_06.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_07.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_08.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_09.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_10.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_11.jpg?raw=true" alt="titlepage" width="800"/>

* Some clarification for the Transform World Matrix: Change the entries in the last column so that the contours are in the center of the viewing volume. **Make sure to use the same transformation for both endo and epi volumes.**  

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_12.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_13.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_14.jpg?raw=true" alt="titlepage" width="800"/>

<img src="https://github.com/mmoth-kurtis/MMotH-Fenics-UK/blob/master/docs/pages/images/FEniCS_LV_segmentation_Page_15.jpg?raw=true" alt="titlepage" width="800"/>

* Run the 1.createLV.py script on the stl files (creating json input for this now to be able to specify output directory and mesh name). This creates .vtk files.
* Run LV_Test.py on the New_Mesh files to yield the HDF5 file. This is what is needed for FEniCS. This step creates the mesh and assigns fiber angles, as well as the local coordinate system for each element.
* The mesh needs the top to be cropped out.
