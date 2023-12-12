# miR_RF APPLICATION
The "miR_RF APPLICATION" repository hosts a machine learning-based application designed to evaluate the authenticity of pre-microRNAs. 

---

# miR_RF APPLICATION Description

The miR_RF application is a predictive application for evaluating pre-microRNAs based on the machine learning algorithm Random Forest. It consists of Python and R scripts designed to process RNAfold Vienna output, extract features, perform machine learning analysis and generate predictions.

### Overview

The miR_RF application is comprised of Python and R scripts:
- **Python Script (Pre-processing):**
  - Extracts features from pre-miRNAs present in RNAfold output files.
  - Converts extracted features into numerical format and produces an intermediate file for R processing.
  - Input: RNAfold output file.
  - Output: Intermediate file with extracted features.
- **R Script (Machine Learning):**
  - Performs machine learning analysis on the intermediate file generated by the Python script.
  - Generates a text file containing predictions for each pre-miRNA in the input.
  - Input: Intermediate file generated by the Python script.
  - Output: Text file with predictions (2 for YES, 1 for NO) for each pre-miRNA.

### Input Requirements

The application accepts RNAfold output files as input, structured in the following format:
- A file with a header line starting with ">" followed by five subsequent lines for each pre-miRNA. 
For every input string, made by the header line followed by the sequence, the output is composed by four lines corresponding to respectively:
1. MFE structure, reported within round brackets;
2. Ensemble structure, with energy reported in square brackets;
3. Centroid structure, with its energy and the minimal base-pair distance to all the structures in
the thermodynamic ensemble, reported in curly brackets.
4. Frequency of MFE structure in ensemble and ensemble diversity
   
- Multi-FASTA format is also supported.

In order to obtain the right input for the miR_RF application, you can install RNAfold Vienna package on your machine, in command line. To install this package run one of the following:

```bash
conda install -c bioconda viennarna
conda install -c "bioconda/label/cf201901" viennarna
```

And then type:

```bash
RNAfold -p -d2 --noLP --noDP --noPS <input_file> > <output_RNAfold_file>
```

This command creates the <output_RNAfold_file> file, which is the input for the miR_RF application. 

- The miR_RF application accommodates a range of input file extensions. Whether it's a .txt, .out, or another format, the application is engineered to process pre-miRNA data 
  effectively, irrespective of the file extension. 
- Important note: the header cannot contain values separated by the tab symbol "\t". Therefore, the application converts by default any "\t" present in the header into a single 
  space " ". 

  For example, from this input:
  
  ```plaintext
  >hsa-let-7a-1  first_example  1
  ```
  to: 
  
  ```plaintext
  >hsa-let-7a-1 first example 1
  ```

  Note: consider that the code only processes and replaces headers presenting the following symbol: "\t". Any other formats will remain in the final output as written in 
  input.  
  
  
### Input Example

Sample input file structure:

```plaintext
>hsa-let-7a-1
UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA
(((((.(((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))))))) (-34.20)
{((((.(((((((((((((((((((((.....(((...((({....}))).))))))))))))))))))))))))))))} [-35.18]
(((((.(((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))))))) {-34.20 d=3.42}
 frequency of mfe structure in ensemble 0.203686; ensemble diversity 5.63
...
```

### Output Example

The output file contains pre-miRNA names and their corresponding predictions:

```plaintext
"miRNA name"       "prediction"
">hsa-let-7a-1"       "2"
```

### Installation

Before beginning the installation, I recommend creating a new directory to neatly store all the requirements for the miR_RF Application. This will facilitate a clearer and more organized environment for running the application efficiently. 
To create a new directory, for example named "miR_RF_application", in your current location, use the following command in the terminal:

```bash
mkdir miR_RF_application
```
This command will create a new directory named "miR_RF_application" within the current location. Users can then put the necessary files here. 

1. Conda Installation in Command Line:
   - Follow the provided instructions in the 'CONDA installation instructions' file to install Conda on your system in the directory just created.

2. Activating the Conda environment:
   - Once Conda is installed, use the provided `configuration_file.yml` file to create an environment suitable for running the miR_RF application.
   Download the `configuration_file.yml` file and copy it in the new directory, as follows:

   ```bash
   cp configuration_file.yml ~/miR_RF_application
   ```
   - In the command line, activate conda with the following command:

   ```bash
   conda activate
   ``` 
   
   - Remain in the directory containing the `configuration_file.yml` file and clone the following command:

   ```bash
   conda env create -f configuration_file.yml
   ```
   Note: Creating the environment may take some time, as Conda downloads and installs the necessary packages and dependencies.

   - After the project environment is created, activate it by using the following command:

   ```bash
   conda activate miR_RF
   ```
   This step ensures that the appropriate environment, complete with all the necessary channels and packages required to run the miR_RF application, is activated. The 
   configuration_file.yml contains a specific set of channel configurations and package installations essential for the execution of the application.
   By following these steps, you will have the correct environment with pre-configured channels and packages, ready to utilize the miR_RF application efficiently.


3. Setting up the directory:
   - Add in the same directory where it is present the `configuration_file.yml` and ".sh" files, the provided following files:
      - `PY_features_extraction.py`: Python script for feature extraction from pre-miRNAs;
      - `make_pred.R`: R script for making predictions using machine learning;
      - `trained_model.RDS`: Includes the pre-trained model data necessary for predictions;
      - `application.py`: Executor program coordinating the feature extraction and prediction processes;
 
   On command line, copy the repository URL in the right directory and write:

   ```bash
   git clone <repository_URL>
   ```
   
   Note: make sure that `configuration_file.yml`, `trained_model.RDS`, `PY_features_extraction.py`, `make_pred.R` and `application.py` are located in the same directory. 
   In order to check, use the ls command along with the file names.
   Type:

   ```bash
   ~/miR_RF_application$ ls
   Anaconda3-2023.09-0-Linux-x86_64.sh configuration_file.yml trained_model.RDS df_feat_ext.py make_pred.R application.py
   ```
  

4. Running the miR_RF application:
   - To utilize the miR_RF application for predicting pre-miRNAs, use the following command in the terminal or command line interface:

   ```bash
   python3 application.py <input_file> <output_file>
   ```

   Replace input_file with the name of the file containing pre-miRNA data in the required format. Similarly, replace output_file with the desired name for the 
   prediction results file.

   Example Usage:

   ```bash
   python3 application.py miRNA_sequences.txt predictions.txt
   ```
   miRNA_sequences.txt: Example input file containing pre-miRNA data.
   
   predictions_output.txt: Output file to store the prediction results.

   Ensure that the input file follows the specified format (see Input requirements). Upon executing this command, the `application.py` program will process the input data, 
   execute feature extraction, and generate predictions using the trained model.


6. Example input files:
   
   Use the provided files, called "miRNA_sequences.txt" and "miRNA_mouse_sequences.txt", if needed, in order to obtain and run an input example.


7. Visualising the output:

   In order to display the final output, you can type "less" or "more" in the command line, for example:

   ```bash
   less <output_file>
   ```
   
   or:

   ```bash
   more <output_file>
   ```
   
