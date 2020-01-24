# MMR
Project for the manual training of the mental rotation test of Jost and Jansen (2019). Contains both code for conducting the experiment as well as code for reading and analyzing the data in R and the resulting dataset.

Code is based on the code at https://github.com/LeonardoJost/MRExperiment

## Conducting the Experiment
The experiment is programmed for use with a steering wheel. The experiment was conducted using a Thrustmaster T150RS and the settings are chosen accordingly.
The experiment contains a pre and post test and three possible training conditions in between and a questionaire after all tests.
### OpenSesame
This folder contains the OpenSesame software (Mathôt, Schreij, & Theeuwes, 2012; https://osdoc.cogsci.nl/) code for the experiment. 
#### Usage
All stimuli and txt files containing instructions and questions should be loaded into the experiment. Parameters can be changed in the inline script 'parameters'. Stimulus properties for each block (models, axes and angles) can be set in 'StimulusPropertiesByBlock'. Joystick settings can be changed in 'InitializeJoystick'. Naming of logfiles is set by ID to avoid manually generating different names and set in 'SetLogfile'. Blocks can be modified in the loop 'Blocks'.

A random stimulus order is generated for each block to account for differing stimulus parameters in each block.

Questions should be entered the following way: (Type)(Options)(ID)(condition). Type should be either: Multiple (for multiple choice), TextInput (for custom text input), ManualCode (requires exact entering of the question ID in a text field to continue), or ShowText (default, only ok button is presented as answer). If type is Multiple, options should include all possible answers separated by ',' (no spaces), for other Types no options should be entered. ID is optional, default is the number of the question. Conditions can be entered if some questions should only be presented depending on previous answers. At the moment only 'ID==answer' is supported and the question is only presented if the question with 'ID' is anwered with 'answer'.

Instructions are shown for a fixed time before allowing continueing to the trials to prevent accidental skipping. The txt may contain a block preceeeded by (Continue), which is shown only after the fixed time and replaced by empty lines (containing single dots as a workaroung for text alignment) before.

The blocks are by default controlled by time and the number of stimuli. If either all stimuli are processed or time is over, the experiment should finish. The transition between blocks is only controlled by time but can be changed to number of stimuli easily in the break if condition of the 'ShowStimuli' loop.

## R Code
R code is in all other files.

R code will read data from files in the specified folder and do some calculations (average participant data, outlier detection for mental rotation data, plotting and averages for mental rotation data) and contains code for basic statistical analysis using linear mixed models (has to be adapted to individual case).