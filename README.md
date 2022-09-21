# Resistance_Evaluation
The ultimate goal of this project is to build a set of interpolated lookup tables, such that each curve in the lookup tables describe the relationship between the on-state resistance and a corresponding temperature for a given drain-current.  Measurements of voltage and current at different junction-temperatures are used to interpolate and extrapolate these relations. 


### Table of Contents

1. [Instructions](#instructions)
2. [Project Motivation](#motivation)
3. [Licensing, Authors, and Acknowledgements](#licensing)

## Instructions: <a name="instructions"></a>

The project was created with Python 3.9.0.
Run the following commands to initiate thw project:

1. create virtual environment in folder **ML_Disaster_Response_Pipeline/**:

  `python3 -m venv env38`

2. activate the virtual environment:

  `source env38/bin/activate`

3. pip install required packages:

  `pip install -r requirements.txt`

4. Use 'Rds_on_eval.py' to produce calculations of resistance ('Rds_Tj_measurement.npz'). 
5. Use 'Bayes Filter.py' to import ('Rds_Tj_measurement.npz') and visualize a sample of the lookup-tables for controllable given drain current.   


## Project Motivation: <a name="motivation"></a>

The optimized lookup tables are used for the accurate and real-time estimation of the junction temperature of a series of silicon-carbide transistor given an accurate measurement of the on-state resistance. 



## Licensing, Authors, Acknowledgements: <a name="licensing"></a>

This corporate message data is from one of the free datasets provided on the
[Figure Eight Platform](https://appen.com/resources/datasets/), licensed under
a [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

Feel free to use my code as you please:

Copyright 2020 Mahmoud Saad

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
