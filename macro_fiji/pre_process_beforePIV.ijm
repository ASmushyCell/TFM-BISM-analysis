// to open this macro, just drag it in FIJI and click on RUN once you have set the correct parameters and selected the proper // active image (the one you want to analyse). 
// macro by Lucas Anger ; 
// this is meant to be used for the TFM

setBatchMode("show");
// note : the analysis is quicker if the data are in 8bit
path_dir = getDirectory("Choose directory where your data are stored");

File.makeDirectory(path_dir);
File.makeDirectory(path_dir + "Input/");
File.makeDirectory(path_dir + "Input/forPIV/");

//concatenate
file = path_dir + "Ref.tif";
if (!File.exists(file)) {
    exit("Error: File 'Ref.tif' not found in the selected folder \n" + "you idiot please try again in the correct location this time");
}
file = path_dir + "WithCells.tif";
if (!File.exists(file)) {
    exit("Error: File 'WithCells.tif' not found in the selected folder \n" + "you idiot + something is wrong like a file is missing ? i can't");
}

open(path_dir + "Ref.tif");
open(path_dir + "WithCells.tif");

selectWindow("WithCells.tif");
InitialStackNumber=nSlices;

run("Stack Sorter");
selectWindow("Ref.tif");
for (i=1; i<InitialStackNumber; i++){
         selectWindow("Ref.tif");
          run("Duplicate...", " ");  
    }
run("Images to Stack", "name=Ref.tif title=[] use");


run("Interleave", "stack_1=Ref.tif stack_2=WithCells.tif");
run("Image Stabilizer", "transformation=Translation maximum_pyramid_levels=1 template_update_coefficient=0.90 maximum_iterations=200 error_tolerance=0.0000001 output_to_a_new_stack");

selectWindow("Ref.tif");
close();
selectWindow("WithCells.tif");
close();
selectWindow("Combined Stacks");
close();
	
//correct luminosity 
selectWindow("Stablized Combined");
run("Illumination correction plugin", "referenceslice=1");
selectWindow("Stablized Combined");
close();
selectWindow("Histogram-matched-Stablized Combined");
run("Image Sequence... ", "dir=" + path_dir + "Input/forPIV/ format=TIFF name=[]");
selectWindow("Histogram-matched-Stablized Combined");
close();

//the end 
showMessage("Process is completed ; clap clap clap this is the end");
