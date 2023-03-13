# Distinguishable Count Algorithm
This program counts the distinguishable structures for one or more RNA
sequences. It detects the symmetry in the order of multiple sequences
and computes the number of distinguishable structures. 


## Requirements

- Boost


## Installation

```bash
git clone git@github.com:masarunakajima/distinguishable_count.git
cd distinguishable_count
make 
```

## Input file format
Input format follows that of [nupack](https://www.nupack.org/). An
example input is shown below.

```
2
GGCTGGTTTCT
GTCTGGGATGCTGGAT
1 2 1 2
```

The first line shows the number `n` of distinct RNA sequences. The
following `n` lines specify the sequences. The last line specifies the
way the sequences are organized. In the above example, the program
counts the distinguishable structures for 4 sequences.
```
GGCTGGTTTCT + GTCTGGGATGCTGGAT + GGCTGGTTTCT + GTCTGGGATGCTGGAT
```

## Executing the program

You can execute the program by running 
```bash
./count_dist examples/input1.in
```
If you want to save the output to a file, you can run
```bash
./count_dist examples/input1.in -o result.out
```
Here is the list of options.

| Option | Description |
| ------ | ----------- |
| -?, -help | Display the help message and exit |
| -o, -out | Specify the output file path |
| -v, -verbose | Enable verbose output |
| -about | print about message |

## Output

The output of 
```bash
./count_dist examples/input1.in
```
would be
```
Number of structures
Total   : 24664534989415
Number of distinguishable structures
Total   : 12332271565072
```

