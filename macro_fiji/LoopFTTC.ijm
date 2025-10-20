//Loop to do FTTC on each PIV files for each frame


//Change parameters here////

args = getArgument();  // Example: "Nframes=n;arg_fttc=pixel=1.32 poisson=0.5 young's=15000 regularization=9e-10 plot=1000 plot=1000 select=[F:\\Stress_defect\\MCF10A\\30kPa_mcf10a_10min_TFM_replicat3\\no13"

// Parse the arguments from matlab 
tokens = split(args, ";");
for (i = 0; i < lengthOf(tokens); i++) {
    pair = split(tokens[i], "~");
    if (pair[0] == "Nframes") Nframes = parseInt(pair[1]);
    if (pair[0] == "arg_fttc") arg_fttc = pair[1];
}


// compute the fttc on a loop for each displacement field 

for (i = 0; i <Nframes; i++){   // 
j=i+1;
run("FTTC ", arg_fttc + "\\Output\\FTTC_Output\\PIV"+j+".txt]");     //***********************************************************
}


