# Manhattan Worlds Odometry (MWO)
This package provides a Matlab implementation of ACCV2016 paper: "Divide and Conquer: Effcient Density-Based Tracking of 3D Sensors in Manhattan Worlds" for the purpose of study.

*   **preprocessing**

    Computing the 3D patches (super point cloud) in advance for efficient experiment.

*   **tracking**

    Core functions of the motion estimation algorithm.
    
*   **evaluation**

    A script for evaluating the result and outputing figures.
    
## Installation

    choose a proper directory and clone by: 
    
    git clone https://github.com/Ethan-Zhou/MWO.git
    
## Usage

*   Download a dataset from https://vision.in.tum.de/data/datasets/rgbd-dataset/download, fr3 cabinet is recommanded.

*   Run preprocessing/Preprocessing.m which gives you 3D patches saved as independent .mat files.

*   Run tracking/main.m which will give you the 3D motion estimation result. 

*   You can tune parameters in load_param_MFVO. Have fun!

## Real-time demos

*  https://www.youtube.com/watch?v=W06ZZki0rTA

## Publication

The approach is descirbed in the following publication:

*  **Divide and Conquer: Effcient Density-Based Tracking of 3D Sensors in Manhattan Worlds** (Yi Zhou, Laurent Kneip, Cristian Rodriguez, Hongdong Li), The 13th Asian Conference on Computer Vision (ACCV 2016), Oral presentation.

You can find the paper at http://users.cecs.anu.edu.au/~u5535909/.

## License

The package is licenced under the MIT License, see http://opensource.org/licenses/MIT.





    
    

    



