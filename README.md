# Model analyzer
This tool allows you to retrieve static information about a given model. You can either supply a file location or a Gurobi model object, and then get in return a dictionary or JSON with the corresponding model information.

This is the code that Kostja wrote 4 years ago, however reformatted by me and with the hope of creating a modular model analysis and productivity tool from it.

You can get it either by cloning this repo or by conda installing it (coming soon); conda is needed since `gurobipy` currently needs conda.
