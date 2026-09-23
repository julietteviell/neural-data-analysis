sript to run in order: 
1. filament_initiation_withdrawal_detection_interactive.m

   ==> from Ramalgo data files (.hf5), it detects the touch of the filament and the withdrawal of the paw automatically.
   
2. classification.m

   ==> processes all single neuron data, generates raw and normalized frequency binning (50ms) around every event.
   ==> classifies neurons by activation or inhibition around every event, for each paw, separated in L/R dorsal horn
   ==> generates heatmaps and neural activity curves in the panel fig4 fig5 and fig6
   
3. correlation_neural_activity_temperature_force.m

   ==> correlates neuronal activity before withdrawal, with the force applied and the thermal gradient.
   ==> !change the first line to "pressure" or "temperature" depending on the test analysis!
   
4. percentage_of_activation_heat_mechanical_ipsi_contra.m

   ==> calculates the percentage of activated cells per bin around an event and generates the percentage of activation figures

data needed:
for the functions: "neuron classification fonctions" must be in the path
"processed_data_classification_script_50msbin.mat" and "m31ramptrials.mat" must be in the path
