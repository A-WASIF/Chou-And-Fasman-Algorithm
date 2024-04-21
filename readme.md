# Chou And Fasman method for secondary structure Prediction of Protien
This code implements the Chou and Fasman method for secondary structure prediction. The method uses amino acid propensity values for alpha helix and beta sheets to predict the secondary structure of a given protein sequence. The code first identifies possible alpha helix and beta sheet sites by scanning the protein sequence. It then filters the possible sites based on the quantity of amino acids with high propensity values for alpha helix and beta sheet formation. The code extends the filtered sites by checking the amino acid propensity values of the neighboring amino acids. The code then assigns the helical and beta sheet regions based on the extended sites. If there are conflicting regions between helix and sheet, the code resolves the conflicts based on the amino acid propensity values. The final resolved secondary structure is printed along with the helical, beta sheet, and conflicting regions.

The code is implemented in Python and can be executed by running the main.py file. The protein sequence is provided as input to the chou_fasman function, which performs the secondary structure prediction and prints the results. The amino acid propensity values for alpha helix and beta sheet are defined within the code. The code also includes a function to print the sequence with the predicted secondary structure as a guide.

To run the code, simply execute the main.py file and provide the protein sequence as input. The predicted secondary structure will be displayed, showing the helical, beta sheet, and resolved regions of the protein sequence.

The Chou and Fasman method is a classic approach to secondary structure prediction and provides valuable insights into the structural characteristics of proteins. The method is based on the observation that certain amino acids have a higher propensity to form alpha helices or beta sheets, which can be used to predict the secondary structure of a protein sequence.

## Terminology Used
Alpha Helix : These are the regions shown by `H` in output which shows that this protien has tendency to belong to alpha helix region.

Beta Strands : These are the regions show by `S` in output which shows that this protien has tendency to belongs to beta strand.

Conflicting Region : These region are where protein shows tendency to belong to alpha helix as well as beta strand, they are shown by `C` in output.

Resolved Region : These are the resolved regions of the conflicting protein regions.

## Output
Predicted region will be shown in following formate :\

Each odd line(starting from 1) will show 50 letters of the given protein sequence and then the just following line (i.e. even ones) will show predicted output.

Symbols used in prediction:

1. Alpha Helix : H
2. Beta Sheet : S
3. Conflicting region : C



## Author
Wasif Ali