The shared code from the paper:"Implementing high performance voltammetry simulation using the implicit parallel algorithm"
http://dx.doi.org/10.1016/j.jelechem.2016.03.030


CONTENTS DESCRIPTION
---------------------------------------------------------------------------------

Folder: 'project'- contains a ready project for modifing the simulation code, requires visual studio 2013 and CUDA 7.5.
   
                  - 'project\hpcv7_5_visualstudio2013.sln' is a visual studio file for lunching the project.
 
                  -  'project\7.5' folder contains the project files including the soruce files:'kernel.cu' and 'kernel_gpu.h'. 
                   
                  -  The classes in the folder 'project\7.5\classes\':'Grid.h', 'Pretomas.h', 'simec.h', 'Specie.h', 'tomas.h', 'World.h'. The classes are needed to run the simulation.
                
                  -  The input file: 'input.exe', and the output file: 'CPU-GPU.txt'.                   
                
                  -  Other files: '7.5.vcxproj', 'cv120.pdb', and the folders 'project\Debug','project\Release' are part of the Visual Studio project's configuration.   

Folder: 'source files' - contains all the source files ('kernel.cu' and 'kernel_gpu.h')
                         and classes ('Grid.h', 'Pretomas.h', 'simec.h', 'Specie.h', 'tomas.h', 'World.h')
                         that are needed to build the programm in different versions of Visual Studio or CUDA.
        
Folder: 'executed.zip' - contains an executeable program for tests in release and debug mode.

--------------------------------------------------------------------------------



INSTRUCTIONS:

*CUDA should be installed in your operating system before running the program.

TEST
-------------------------------------------------------------------------------

- To test the program, there is no need for Visual Studio, but CUDA should be installed in your computer.

- The executed files can be found in the folder 'executed' in the release mode and in the debug mode (If you do not have Visual Studio installed, use the release mode).

- To run the program in the release or debug mode, click on the '7.5.exe' file in the sub-folder ...\executed\Release\7.5.exe or ...\executed\Debug\7.5.exe.

- The program uses the input parameters from the file 'input.txt', which is in the same folder of the executing file '7_5.exe'. Edit this file in order to change the simulation parameters.

- The executed program runs the same simulation on the CPU and following that on the GPU. 
  An output file with the name 'CPU-GPU.txt' is created when the program is finished.
  The  file is located on the same folder of the executable file '7.5.exe', and contains four columns: The time, voltage, GPU flux, CPU flux. 

------------------------------------------------------------------------------


BUILD AND MODIFY THE CODE
-----------------------------------------------------------------------------
###  The project was programmed in C++ with CUDA, and a Visual Studio 2013  project have been set.
   The ready project  is located on in the folder 'project' and can be lunched by clicking on the file 'hpcv7_5_visualstudio2013.sln'.  
   In order to be able to run this project, Visual Studio 2013 and CUDA 7.5 should be installed. 

###  In the case that a later version of CUDA is installed or other version of Visual Studio (should support CUDA), a new project should be created as follows:
   - Open Visual Studio, and navigate to File->New -> Project.
   - Look for NVIDIA tab in the project type, and choose your installed CUDA runtime environment. 
   - Create the project in a new folder.
   - Copy the folder 'source file' to the new folder you have created for the project.
   - In Visual Studio, choose add->Existing Item... 
   - Navigate inside the copied 'source file' folder and add the source files kernel.cu, and kernel_gpu.h.
   - Build/Debug the solution.

###  In a different IDE or different OS the user can find all the source files in the folder 'source files'.

----------------------------------------------------------------------------





  
