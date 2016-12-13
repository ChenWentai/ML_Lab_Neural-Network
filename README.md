# Machine-Learning-Lab-Neural-Network
Intelligent Systems - Fall 2016
2nd Lab Session's Report
By Wentai Chen  and Di Huang
Part 1: Design a neural network<br>
=========
Describe how you have computed the network architecture and weights. 
You should describe the network according to the following format: 

Input layer:  2 units, x1   x2
First hidden layer:
    2 neurons:
	Neuron 1: w11 = <2>
                  w21 = <3>
                  b1 = <230>
	Neuron 2: w12 = <1.1>
                  w22 = <-1>
                  b2 = <0>
	
Output layer:
    3 neurons:
	Neuron 1: w11 = <1>
                  w21 = <1>
                  b1 = <1.5>
	Neuron 2: w12 = <1>
                  w22 = <-1>
                  b2 = <1.5>
	Neuron 3: w13 = <0>
                  w23 = <-1>
                  b3 = <-0.5>
Of course, it is your task to define the number of layers, the number of neurons per layer, and the exact values for the weights. 
Part 2: Implementation of the MLP simulator<br>
=======
Provide source code and explanations (comments) for both the feedforward() and backpropagation() procedure. 
The feedforward function.
The backpropagation function.
Part 3: Training and Recall experiments<br>
=======
Detail and comment all the experiments achieved: 
•	Network training, 
•	Recall performance (on the training set set30.x1x2rgb and on the recall set set120.x1x2rgb 
•	The effect of the learning rate, the number of training cycles and the RMS error stopping criterion. Here you could include a graph showing the RMS error as a function of the number of training cycles for different learning rate values. 
