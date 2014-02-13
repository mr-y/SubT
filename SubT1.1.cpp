/********************************************************************
Copyright (C) 2010 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: kryberg@utk.edu
*********************************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <algorithm>

using namespace std;

/*** Function declarations ***/
void mcmc(double *nodeheight, int ntaxa, int nsubstitute, int ngen, int sampfreq, char yule, double priorNetdiv, int priorShape, double window, double windowNode, int treeno, char *nodelog);
double bdLh(double *nodeheight, double b, double d, int ntaxa);
void MCMC_summarize (istream& infile, int burnin);
void help();

/*** Probability distributions ***/
double gamma(double x, int k, double b);
double uniform(double x, double a, double b);
double normal(double a);
double normalprob (double x, double mean, double sd);
double exponential (double lamda);

/*** Math ***/
int factorial(int x);

/*** Classes ***/

/*** Class for phylogenetic trees ***/
class tree {
    public:
        struct node { //store the information for each node of the tree
            double branchlength; //branchlength of branch leading to node
            string nodelable; //name of node
            node *left; //left daughter node
            node *right; //right daughter node
            node *parent; //the parent of the node
        };
        tree () { root = new node; }; //initilize the clase by creating the root
        ~tree  () { destroy_tree (root); }; //delete the class starting from the node
        //function to get node heights, destroys the tree
        double *node_heights ( double nodeheights[], int n_nodes ); 
        //parses a newik formated file and stores it in the tree (object)
        void tree_file_parser( istream& infile );

    private:
        tree (const tree&, tree tree_copy); //can not copy tree objects (yet)
        //deletes a part of the tree from given node
        void destroy_tree (node *leaf);
        node *root; //stores the location of the root
};

/*** Class to handel data as a binary tree ***/
class data_tree {
    public:
        struct node { //store the information
            double value; //store a floating point value
            node *left; //points to lower values
            node *right; //point to higher values
        };
        data_tree () { root = 0; total = 0; number_of_values = 0; } //initilize the values
        ~data_tree () { destroy_tree(root); } //deletes the tree starting from root
        void insert_data (double x); //function ti insert data
        double mean() { return total/number_of_values; } //calculate mean
        double median() { //get the median value
            if ( n_nodes() % 2 == 0 ) return (get_value((n_nodes()/2)) + get_value((n_nodes()/2)+1))/2; //if even number of nodes get the mean of the equaly middle values
            else return get_value((n_nodes()/2)+1); //get middle vale
        }
        double sum() { //sum up values in tree
            if (root == 0) return 0; //if the root is empty there is no sum
            return sum(root); //call sum function
        }
        int n_nodes() { //get number of nodes
            if (root == 0) return 0; //if the root is empty there is no values to sum up
            return n_descendants(root); //calculate the sum from the root up
        }
        double sd() { //calculate standard deviation
            if ( root == 0) return 0; //if root is empty there is no values
            return sqrt(sum_square(root)/n_nodes()); //take the square root of the mean of sum of the square differences between values and the mean
        }
        double max() { //get the highest value in the tree
            node *present = root; //start at root
            if (root == 0) return 0; //if root emty there are no values
            while (present->right != 0) present=present->right; //while not at a tip go to the right node
            return present->value; //the highest value is in the right most node
        }
        double min() { //get the lowest value in the tree
            if ( root == 0 ) return 0; //if the root is empty there are no values
            node *present = root; //start at root
            while (present->left != 0) present=present->left; //while not at a tip go to the left
            return present->value; //the left most tip has the lowest value
        }
        double get_value( int i ); //function to get a value counted from the lowest (value 1 beeing the lowest)

    private:
        data_tree (const data_tree&, data_tree vector_copy); //objects of the class can not be copied
        void destroy_tree ( node *leaf ) { //delete the tree from the given node
            if(leaf!=0) { //while not at a tip delete next nodes
                destroy_tree ( leaf->left ); //delete node to the left
                destroy_tree ( leaf->right ); //delete node to the right
                delete leaf; //delete present node
            }
        }
        void insert_data ( double x, node *leaf ); //inser data somwhere after given node
        double sum(node *leaf) { //sum the values in the tree from the given node
            if ( leaf == 0 ) return 0; //if leaf is empty return 0
            return leaf->value + sum(leaf->left) + sum(leaf->right); //return the sum of the value of the leaf and all values to the left and right
        }
        double sum_square ( node *leaf ) { //calculate the sum of the square of the difference between the mean and the value for all nodes higher in the tree
            if ( leaf == 0 ) return 0; //if the leaf is empty return 0
            return (( leaf->value-mean() ) * ( leaf->value-mean() )) + sum_square(leaf->left) + sum_square(leaf->right); //return the sum of the square difference of this leaf and all nodes to the left and right
        }
        int n_descendants ( node *leaf ) { //returns the number of daughter nodes for a node +1
            if (leaf == 0) return 0; //if a empty leaf return 0
            return 1 + n_descendants(leaf->left) + n_descendants(leaf->right); //return 1 plus the number of descendants to the left and right
        }
            
        node *root; //store the root location of the tree
        double total; //store the total off all values in the tree, needed to make mean() fast so sd() is resonable fast
        int number_of_values; //store the total number of nodes in the tree for the same reason total is stored
};

/** Functions for sorted vector that ***
*** should not be inlined            **/
void data_tree::insert_data ( double x ) { //function to insert data
    if ( root == 0 ) { //if the node has not been initiated do so
        root = new node; //allocate memory for the root
        root->value = x; //store the value
        root->left = 0; //set the left leaf to a null pointer
        root->right = 0; //set the right leaf to a null pointer
    }
    else { insert_data( x, root ); } //otherwise insert data, starting to look for where at the root

}

void data_tree::insert_data ( double x, node *leaf ) { //insert data somwhere after the given node
    if ( x < leaf->value ) { //if the value is less than the value in the present leaf it belonges to the left
        if ( leaf->left != 0 ) insert_data ( x, leaf->left); //if not at a tip proceed (recursive function)
        else { //if at tip create node
            leaf->left = new node; //allocate memory for the left node
            leaf->left->value = x; //insert value
            leaf->left->left = 0; //set left leaf to null pointer
            leaf->left->right = 0; //set right leaf to null pointer
        }
    }
    else if ( x >= leaf->value ) { //if value higher (or equal) to present node insert it to the right
        if ( leaf->right !=0 ) insert_data ( x, leaf->right ); //if not at a tip proceed (recursive finction)
        else { //if at tip insert value
            leaf->right = new node; //allocate memory for right node
            leaf->right->value = x; //insert value
            leaf->right->left = 0; //set left leaf to null pointer
            leaf->right->right = 0; //set right leaf to null pointer
        }
    }
    total += x; //add value to the total
    ++number_of_values; //add one to the number of nodes
}

double data_tree::get_value( int i ) { //get the i:th lowest value
    if ( root == 0 ) { cerr << "Trying to get value from uninitiated data_tree; returning 0." << endl; return 0; } //if empty tree return 0, this is not strictly correct as it should return nothing
    if ( i < 0 || i > n_nodes() ) { cerr << "get_value function in data_tree (class) called with value out of range; returning 0." << endl; return 0; } //if the node is out of range return 0, should return nothing
    node *present = root; //start at root
    while ( i > 0 ) { //loop, should stop when reaching right value but if i is < 0 something is wrong
        if (i <= n_descendants(present->left) ) present = present->left; //if there is more nodes than i nodes to the left the value we are looking for is to the left
        else if ( i == n_descendants(present->left) + 1 ) return present->value; //when we have i-1 nodes to the left we are at the value we are looking for
        else if ( i > n_descendants(present->left) ) { //if i is higher than the number of nodes to the left we are looking for a higher value
            i -= n_descendants(present->left) + 1; //after we moved to the right we will have excluded the number of nodes to the left plus one from our search (as lower values)
            present = present->right; //the value we are looking for is higher than the present node and nodes to the left so ove to the right
        }
    }
    return present->value; //should not get here so shuld probably give error message and return 0 instead, but the value we are at should be closer to the real value we are looking for
}

/***********************************
************************************
*** Main function of the program ***
************************************
***********************************/
main (int argc, char *argv []) {
   char *filename;          //Pointer to stored file name
   char temp[100];          //Temporary variable to hold imput
   char yule='Y';           //Keep track of which model to use
   int ntaxa=0;             //Variable to store number of genera
   int i=0;                 //general counter
   int treeno=0;            //To keep track of what tree is being analyzed
   int nsubstitute=0;       //Number of taxa not included phylogeny and whose node heights should be estimated
   double window=0.1;       //Window width for proposition of new state of diversification rates in mcmc
   double windowNode=0;     //Window width for proposition of new state of nede heights in mcmc
   double priorNetdiv=0.1;  //Prior mean for the net diversification rate
   int priorShape=3;        //Shape of gamma distribution for net diversification
   int ngen=1000000;        //Variable for number of generations to run
   int sampfreq=1000;       //Variable for sample frequency
   char nodelog[100]="NOLOGFILE"; //string to keep log file name for node heighs
   int bracket_balance=0;   //Variable to keep count on the balance of right and left square brackets
   srand(time(NULL));       //Set the random seed by clock

   //Print the cammand used to start the program
   cout<< "The program was called using the following command: ";
   while (argv[i]) cout << argv[i++] << " ";
   cout << endl << endl;

   //See to that there are the right amount of arguments
   if (argc>2) {
      int burnin = 0;
      for (i=1; i < argc; ++i) {
         if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--summerize")) { filename = argv[++i]; }
         else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--burnin")) { burnin = atoi(argv[++i]); }
         else {cout << "Do not recognize arguments, just give name of the infile, -s and name of MCMC file, or --help (-h).\n"; return 1; }
      }
      ifstream infile (filename);
      // If the file will not open quit:
      if (!infile.is_open()) { cout << "Could not open \"" << filename << "\". Quitting.\n"; return 1;}
      MCMC_summarize(infile, burnin); //start summarizing
      infile.close();
      return 0;
   }

   if (argc==1) {cout << "Too few arguments, give name of the infile or --help (-h).\n"; return 1;}

   //If asked for print help
   if (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")) { help(); return 1; }

   //Read infile name
   else if (argc==2) { filename = argv[1];} //if only one argument, beside program name, assume it's a file name

   ifstream infile (filename);
   // If the file will not open quit:
   if (!infile.is_open()) { cout << "Could not open \"" << argv[1] << "\". Quitting.\n"; return 1;}

   // If open pars the infile:
   while (infile) {
       // Read next word:
       infile >> temp;
       /*** Remove comments in square brackets ***/
       i=0; //Set counter to 0
       while ( temp[i] != '\0') { //read each character till end of string
          if (bracket_balance < 0) { cout << "Found a unmatched ]. Square brackets need to be balanced. Quitting." << endl; return 1; }
          else if ( temp[i] == '[' ) { //Start of comment to be removed
             bracket_balance++;   //Note one [ for balance
             int j=i;
             while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
             continue; //check next character
          }
          else if ( temp[i] == ']' ) { //End of comment
             bracket_balance--;        //Note one ] for balance
             int j=i;
             while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
             continue; //check next character
          }
          else if (bracket_balance > 0) { //Means we are still in a comment
             int j=i;
             while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
             continue; //check next character
          }
          i++; //If not square bracket or between square brackets check next character
       }
       /******************************************/
       if (temp[0]=='\0') continue; //If nothing in string read next word

       else if (!strcmp(temp,"start") || !strcmp(temp,"start_tree") || !strcmp(temp,"start_file")) {break;}

       // Read number of taxa:
       else if (!strcmp(temp,"ntaxa")) {
         infile >> temp;
         ntaxa=atoi(temp);
         }
       // Read number of substitute taxa:
       else if (!strcmp(temp,"nsubstitute")) {
         infile >> temp;
         nsubstitute=atoi(temp);
       }
       // Read modifier for how wide the suggestion window is in the MCMC for diversification parameters:
       else if (!strcmp(temp,"windowidth") || !strcmp(temp,"windowwidth") || !strcmp(temp,"window")) {
         infile >> temp;
         window=atof(temp);
       }
       // Read modifier for how wide the suggestion window is in the MCMC for node heights:
       else if (!strcmp(temp,"nodewindow")) {
         infile >> temp;
         windowNode=atof(temp);
       }
       // Read mean for net diversification prior distribution:
       else if (!strcmp(temp,"meanGamma")) {
         infile >> temp;
         priorNetdiv=atof(temp);
       }
       else if (!strcmp(temp,"shapeGamma")) {
         infile >> temp;
         priorShape=atoi(temp);
       }
       // Read number of generations to run the MCMC:
       else if (!strcmp(temp,"ngen")) {
         infile >> temp;
         ngen=atoi(temp);
       }
       // Read how often to sample in the MCMC:
       else if (!strcmp(temp,"sampfreq")) {
         infile >> temp;
         sampfreq=atoi(temp);
       }
       //Read what model to use Yule (Y) or birth-death (anything else)
       else if (!strcmp(temp,"model")) {
         infile >> temp;
         yule=temp[0];
       }
       //Read in file name for logfile for node heights of substitute taxa
       else if (!strcmp(temp,"nodelog")) {
         infile >> temp;
         i=0;
         do {nodelog[i]=temp[i++];}
         while (temp[i] != '\0' && temp[i] != '\n' && temp[i] != '\t' && temp[i] != ' ' && i<100);
         nodelog[i]='\0';
       }
       else {
          cout << "Could not understand \"" << temp << "\" in the input file. Please check the manual for available commands." << endl;
          return 1;
       }
   }
   // Close infile again:
   infile.close();
   // If infile still open give warning and quit:
   if (infile.is_open()) { cout << "Could not close file!!! Quitting.\n"; return 1;}

   //If the square brackets are unbalanced print warning and quit
   if (bracket_balance > 0) { cout << bracket_balance << " too many [. Square brackets need to be balanced. Quitting." << endl; return 1; }
   else if (bracket_balance < 0) { cout << bracket_balance * -1 << " too many ]. Square brackets need to be balanced. Quitting." << endl; return 1; }

   // Set windowNode if not set and there are substitute taxa
   if (nsubstitute !=0 && windowNode <= 0) {windowNode = window;}

   // Compensate for lowercase letter when giving yule model
   if (yule=='y') yule='Y';
   if (yule=='d' || yule=='b' || yule=='D') yule='B';

   // Test for variables:
   if (ntaxa<1) { cout << "Did not find number of taxa (ntaxa). Number of taxa must be at least 1 and given in the beginning of the file by \"ntaxa\" followed by the number of taxa. Quitting.\n"; return 1;}
   if (ngen<1) { cout << "Number of generations (ngen) < 1. Must be at least one. Quitting." << endl; return 1;}
   if (sampfreq<1 || sampfreq>=ngen) { cout << "Sample frequency (sampfreq) must be larger than zero and smaller than the number of generations. Quitting." << endl; return 1;}
   if (yule!='Y' && yule!='B') { cout << "The model (model) must be set to y, Y, yule, or Yule (for Yule [pure birth mode] model), or d, birth-death, or Birth-Death. Quitting." << endl; return 1;}
   if (priorNetdiv<=0) { cout << "The net diversification prior must have a mean (meanGamma) larger than zero. Quitting." << endl; return 1;}
   if (priorShape<=0) { cout << "The shape of the net diversification prior (shapeGamma) need to be larger than zero. Quitting." << endl; return 1;}
   if (window<=0) { cout << "The window width multiplier (windowidth) must be larger than zero. Quitting." << endl; return 1;}
   if (nsubstitute < 0) { cout << "The number of substitute taxa (nsubstitute) can not be a negative number. Quitting." << endl; return 1;}
   if (nsubstitute > 0 && strcmp(nodelog,"NOLOGFILE")) {
      ofstream nodeLog;
      nodeLog.open (nodelog);
      if (!nodeLog.is_open()) { cout << "Could not open node log file (nodelog) \"" << nodelog << "\". Quitting.\n"; return 1;}
      nodeLog.close();
      if (nodeLog.is_open()) { cout << "Could not close node log file (nodelog) \"" << nodelog << "\". Quitting.\n"; return 1;}
   }

   //Print the settings for the MCMC
   cout << "Number of species there are node heights for (ntaxa): " << ntaxa << endl;
   cout << "Number of substitute taxa (nsubstitute): " << nsubstitute << endl;
   cout << "Number of generations to run the MCMC (ngen): " << ngen << endl;
   cout << "Frequency with which the MCMC is sampled (sampfreq): " << sampfreq << endl;
   cout << "The model used (model): ";
   if (yule=='Y') cout << "Yule" << endl;
   else cout << "birth-death" << endl;
   cout << "The mean of the prior distribution for the net diversification rate (meanGamma): " << priorNetdiv << endl;
   cout << "The shape of the gamma distribution used as prior for the net diversification rate is (shapeGamma): " << priorShape << endl;
   cout << "The window width multiplier for net diversification rate and relative extinction rate (windowidth): " << window << endl;
   cout << "The window width multiplier for substitute taxa node height (nodewindow): " << windowNode << endl;
   if (strcmp(nodelog,"NOLOGFILE")) {
      cout << "Node heights for substitute taxa are written to (nodelog): ";
      int k=0;
      while (nodelog[k]!='\0') {cout << nodelog[k++];}
   }

   // Create an array to store the node ages:
   double *nodeages= new double[ntaxa-1];
   // Open infile or quit:
   infile.open (filename);
   if (!infile.is_open()) { cout << "Could not open \"" << argv[1] << "\". Quitting.\n"; return 1;}

   // Sett counter
   i=-1; //counter <0 mean not ready to read data
   bracket_balance=0;
   char data_as_trees = 'n'; //flag to see if trees or node depths will be read
   // Read infile tree by tree and execute MCMC
   while (infile) {
      infile >> temp;
      /*** Remove comments in square brackets ***/
      int k=0; //Set counter to 0
      while ( temp[k] != '\0') { //read each character till end of string
         if (bracket_balance < 0) { cout << "Found a unmatched ]. Square brackets need to be balanced. Quitting." << endl; return 1; }
         else if ( temp[k] == '[' ) { //Start of comment to be removed
            bracket_balance++;   //Note one [ for balance
            int j=k;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         else if ( temp[k] == ']' ) { //End of comment
            bracket_balance--;        //Note one ] for balance
            int j=k;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         else if (bracket_balance > 0) { //Means we are still in a comment
            int j=k;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         k++; //If not square bracket or between square brackets check next character
      }
      /******************************************/
      if (temp[0]=='\0') continue; //If nothing in string read next word

      /*** this option will prepare to read nodeheight data ***/
      if (!strcmp(temp,"start")) { // time to read data
         i=0; // set counter to 0 means ready to read data
         treeno=1; //start with tree one (per dfinition)
         cout << endl << endl;
         cout << "Starting MCMC" << endl << "*******************************************************" << endl << endl;
      }
      /*** This option will read tree data within the input file ***/
      if (!strcmp(temp,"start_tree")) { 
         i=0; // set counter to 0 means ready to read data
         treeno=1; //start with tree one (per dfinition)
         data_as_trees = 'y';
         cout << endl << endl;
         cout << "Starting MCMC" << endl << "*******************************************************" << endl << endl;
         while (1) {
             tree intree;
             for (int j=0; j < ntaxa-1; ++j) { nodeages[j] = 0; }
             intree.tree_file_parser(infile);
             nodeages = intree.node_heights(nodeages,ntaxa-1);
             if(nodeages[0] < 0.00000000000001) {
                 if (treeno == 1) {
                     cout << "No trees read!!! This may be due to error in tree format (including no branch length or zero tip branch length) or miss specification of ntaxa." << endl;
                 }
                 break;
             } 
             mcmc(nodeages, ntaxa, nsubstitute, ngen, sampfreq, yule, priorNetdiv, priorShape, window, windowNode, treeno, nodelog);
             treeno++; //Done with present tree so next tree is in line            
         }
      }
      /*** This option will open a file and read trees from that file ***/
      if (!strcmp(temp,"start_file")) {
         i=0; // set counter to 0 means ready to read data
         treeno=1; //start with tree one (per dfinition)
         data_as_trees = 'y';
         infile >> temp; //read treefile name
         ifstream treefile (temp); //open tree file
         // If the file will not open quit:
         if (!treefile.is_open()) { cout << "Could not open treefile \"" << temp << "\". Quitting.\n"; return 1;}         
         else { cout << endl << "Reading trees from \"" << temp << "\"."; }
         cout << endl << endl;
         cout << "Starting MCMC" << endl << "*******************************************************" << endl << endl;
         while (1) {
             tree intree;
             for (int j=0; j < ntaxa-1; ++j) { nodeages[j] = 0; }
             intree.tree_file_parser(treefile);
             nodeages = intree.node_heights(nodeages,ntaxa-1);
             if(nodeages[0] < 0.00000000000001) {
                 if (treeno == 1) {
                     cout << "No trees read!!! This may be due to error in tree format (including no branch length or zero tip branch length) or miss specification of ntaxa." << endl;
                 }
                 break;
             }    
             mcmc(nodeages, ntaxa, nsubstitute, ngen, sampfreq, yule, priorNetdiv, priorShape, window, windowNode, treeno, nodelog);
             treeno++; //Done with present tree so next tree is in line            
         }
      treefile.close();
      }
      if (!strcmp(temp,"end") || !strcmp(temp,"end;") || !strcmp(temp,"End") || !strcmp(temp,"End;") || !strcmp(temp,"END") || !strcmp(temp,"END;")) { //end reading data and executing MCMC
         if (i>0 && i< ntaxa-1 && data_as_trees == 'n') cout << "Number of node heights less than ntaxa -1. Missing data. Last set of node heights were not analyzed." << endl; //have not read all data in last dataset as defined by ntaxa
         break; 
      }
      if (atof(temp)>0.0 && i>-1 && data_as_trees == 'n') {
         nodeages[i++] = atof(temp); 
         if (nodeages[i-1] <= 0.0) { cout << "All node ages must be > 0. Node age " << i-1 << " of tree no " << treeno << " was found to be: " << nodeages[i-1] << endl; return 1;} // does not seem to read negative values but still testing for them
         if (i == ntaxa-1) { // start MCMC if a full set of data has been read
            i=0; 
            mcmc(nodeages, ntaxa, nsubstitute, ngen, sampfreq, yule, priorNetdiv, priorShape, window, windowNode, treeno, nodelog);
            treeno++; //Done with present tree so next tree is in line
         }
      }
   }
   infile.close(); // Close infile
   cout << endl; // Make sure the comand prompt get a new line

/*** Clean up ***/
delete [] nodeages;

}


void mcmc (double *nodeheight, int ntaxa, int nsubstitute, int ngen, int sampfreq, char yule, double priorNetdiv, int priorShape, double window, double windowNode, int treeno, char *nodelog) {
   double Lh=0; //The Likelihood of the parameters
   double PriorDiv=0; //The prior of the diversification parameters
   double PriorNode=0; //The prior of the node heights
   double propLh; //The Likelihood of the proposed parameters
   double propPriorDiv=0; //The prior of the proposed Diversification parameters
   double propPriorNode=0; //The prior of the proposed node heights
   int nspecies=ntaxa+nsubstitute; //Total number of species in the group including substitute taxa
   double netdiv=log(nspecies)/nodeheight[0]; //The net diversification rate
   double extprop=0.1; //The extinction proportion
   double propnetdiv; //The proposed net diversification rate
   double propextprop=0; //The proposed extinction proportion
   double nodeandsubstitute[ntaxa+nsubstitute-1]; //An aaray to store the node heights of all species
   double propnodeandsubstitute[ntaxa+nsubstitute-1]; //An aaray to store the proposednode heights of all species
   double substitute[nsubstitute]; //An array to store the node heights of the substitute taxa 
   double propsubstitute[nsubstitute]; //An array to store the proposed node heights of the substitute taxa
   int acceptence=0; //To store the number of accepted proposals
   int acceptenceNode=0;
   int x=0; //counter
   int i=0; //counter
   double rest=0; //to keep the difference between node heights
   //Open file for node heights
   ofstream nodeLog;
   if (nsubstitute != 0 && strcmp(nodelog,"NOLOGFILE")) {
      nodeLog.open (nodelog,ios::app);
   }
   if (yule == 'Y') { extprop=0; }

/*** Define reasonable starting node heights for substitute lineages ***/ 
   nodeandsubstitute[nspecies-1]=1/double(nspecies);
   for (i=nspecies-2; i>=0; i--) { nodeandsubstitute[i] = nodeandsubstitute[i+1]+(1/(double(i)+1)); }
   for (i=nspecies-1; i>=0; i--) { nodeandsubstitute[i]/=nodeandsubstitute[0]; }
   // Remove the suggested node heights closest to the observed node heights
   for (i=0; i<ntaxa-1; i++) { 
      rest=nodeheight[0]+1;
      for (int j=0; j<nspecies-1; j++) {
         if (fabs((nodeandsubstitute[j]*nodeheight[0])-nodeheight[i])<rest) {
            rest=fabs(nodeandsubstitute[j]*nodeheight[0]-nodeheight[i]); 
            x=j;
         }
      }
      nodeandsubstitute[x]=10; // Give a unrealisticlly high value closest to the observed value
   }
   x=0; // reuse x as an counter
   // Keep suggested substitution taxa with realistic values
   for (i=0; i<nspecies-1; i++) { if (x<nsubstitute && nodeandsubstitute[i]<=1) { substitute[x++]=nodeandsubstitute[i]; }} 
/**********************************************************************/

   for (i=0; i<nspecies-1; i++) { nodeandsubstitute[i]= 0; } // Set nodeheights to 0
   for (i=0; i<ntaxa-1; i++) { nodeandsubstitute[i]=nodeheight[i];  propnodeandsubstitute[i]=nodeheight[i];} // Input observed node heights
   for (i=0; i<nsubstitute; i++) { nodeandsubstitute[i+ntaxa-1]= substitute[i]*nodeheight[0]; } // Input substitute node heights

   // If first tree print headers for the Markov chain, otherwise don't:
   if (treeno == 1) {
      cout << "Tree_no\t" << "Generation\t" << "Log_Lh\t" << "Prior_P\t" << "Post_P\t" << "Netdiv\t" << "Rel_ext_rate\t" << "Acceptance_rate";
      if (nsubstitute != 0) { cout << "\tAcceptance_rate_nodes";}
      cout << endl;
   }
   Lh = bdLh(nodeandsubstitute, netdiv, extprop, nspecies); // Calculate log likelihood
   PriorDiv = log(gamma(netdiv, priorShape, priorNetdiv) * uniform(extprop,0,1)); // Calculate log prior probability
   for (int j=0; j < nsubstitute; j++) { PriorNode += log(uniform(substitute[j],0,1));  }


   // Print start values:
   for (i=0; i < nsubstitute; i++) { PriorNode = PriorNode + log(uniform(substitute[i],0,1));  }
   cout << treeno << "\t0\t" << Lh << "\t" << PriorDiv+PriorNode << "\t" << Lh+PriorDiv+PriorNode << "\t" << netdiv << "\t" << extprop << "\t" << acceptence;
   if (nsubstitute != 0) { cout << "\t" << acceptenceNode;}
   cout << endl;

/*******************/
/*** Start chain ***/
/*******************/
   for (i=1; i<=ngen; i++) {

/**** First suggest new diversification rates and test if to accept them ****/
/*** Suggest new diversification rates ***/
      propnetdiv=netdiv-normal(window); // Suggest new net diversification rate
      if (yule != 'Y') {propextprop=extprop-normal(window);} // If not Yule model uppdate extinction proportion

      propLh=bdLh(nodeandsubstitute, propnetdiv, propextprop, nspecies); // Calculate new log likelihood
      propPriorDiv= log(gamma(propnetdiv, priorShape, priorNetdiv) * uniform(propextprop,0,1)); // Calculate new log posterior probability


      //Test if accepting
      if ( exp((propLh+propPriorDiv)-(Lh+PriorDiv)) > (rand()/(double)RAND_MAX)) {
         Lh=propLh;
         PriorDiv=propPriorDiv;
         netdiv=propnetdiv;
         extprop=propextprop;
         acceptence++;
      }
/************************************************************/

/**** Propose new node heights and test if to accept them ****/
      // Uppdate selected node height
      if (nsubstitute != 0) {
         double proposalratio=1; 
         for (int j=0; j < nsubstitute; j++) { propsubstitute[j] = substitute[j]-normal(windowNode/double(nspecies)/*substitute[j]*/);}
         for (int j=0; j < nsubstitute; j++) { propnodeandsubstitute[j+ntaxa-1] = propsubstitute[j]*nodeheight[0]; }
      
         propLh=bdLh(propnodeandsubstitute, propnetdiv, propextprop, nspecies); // Calculate log likelihood for proposed node heights
         // Calculate new log posterior probability
         propPriorNode=0;
         for (int j=0; j < nsubstitute; j++) { propPriorNode += log(uniform(propsubstitute[j],0,1));  }

         // To accept or not accept proposed state? If accept uppdate parapeters and propbabilities

         // Acceptence criteria if having substitute taxa
//         for (int j=0; j < nsubstitute; j++) {proposalratio*=(normalprob((propsubstitute[j]-substitute[j])/(window*propsubstitute[j]),0,1)/normalprob((propsubstitute[j]-substitute[j])/(window*substitute[j]),0,1));}
         if ( exp((propLh+propPriorNode)-(Lh+PriorNode))/*proposalratio*/ > (rand()/(double)RAND_MAX)) {
            Lh=propLh;
            PriorNode=propPriorNode;
            for (int j=0; j < nsubstitute; j++) { substitute[j]=propsubstitute[j]; nodeandsubstitute[j+ntaxa-1] = propnodeandsubstitute[j+ntaxa-1]; }
            acceptenceNode++;
         }
      }
/************************************************************/

      // If time to sample the chain do so:
      if (fmod(i,sampfreq)==0) {
         cout << treeno << "\t" << i << "\t" << Lh << "\t" << PriorDiv+PriorNode << "\t" << Lh+PriorDiv+PriorNode << "\t" << netdiv << "\t" << extprop << "\t";
         if (nsubstitute != 0) {
            cout << (double)acceptence/i << "\t" << (double)acceptenceNode/i << endl;
            nodeLog << "[treeno:_" << treeno << ",_generation:_" << i << "] ";
            for (int j=0; j < nsubstitute; j++) { nodeLog << nodeandsubstitute[j+ntaxa-1] << " "; }
            nodeLog << endl;
         }
         else {cout << (double)acceptence/i << endl;}
      }

   }
/******/
//Close the logfile for node heights
if (nsubstitute != 0) { nodeLog.close(); }

}
 
/*** Returns log likelihood for bd model ***/
double bdLh(double *nodeheight, double netdiv, double extprop, int ntaxa) {
   double lik=0;
   int i;
   double sum=0;
   double product=log(1/pow(exp(netdiv*nodeheight[0])-extprop,2));
   for (i=1; i<ntaxa; i++) { 
      sum += nodeheight[i];
      product += log(1/pow(exp(netdiv*nodeheight[i])-extprop,2));
   }
   lik=(netdiv*sum)+ntaxa*log(1-extprop)+product;
   for (i=ntaxa-1; i>1; i--) { lik=lik+log(i*netdiv); }
   return lik;

}

/*** Distributions ***/
/* Probability distributions */
double gamma(double x, int k, double b) { //Takes value for which to calculate probability (x), shape (k), and inv. scale (b).

   return pow(x,k-1)*(exp(-x/(1/b))/(pow(1/b,k)*factorial(k-1)));

}

double uniform (double x, double a, double b) {
   if (a > b) return 0;
   if (x < a) return 0;
   if (x > b) return 0;
   else return 1/(b-a);
}

double normalprob (double x, double mean, double sd) {
   
   return (1/sqrt(2*3.14159265*pow(sd,2)))*exp(-(pow(x-mean,2)/(2*pow(sd,2))));

}

/***/

/* Random distributions */
//Generate normal distributed random value using the Marsaglia polar method
double normal(double a) {
   double s=2;
   double x;
   double y;

   while (!(s<1)) {
      x = (rand()/((double)RAND_MAX/2))-1;
      y = (rand()/((double)RAND_MAX/2))-1;
      s = (pow(x,2)+pow(y,2));
   }
   return x*sqrt(-2*log(s)/s)*a;

}
/***/
double exponential (double lambda) {

   return -log(rand()/(double)RAND_MAX)/lambda;

}
/***/
/*** Math ***/
int factorial(int x){
 int sign=1; // to correct for if the number is posetive or negative
 int i;

 if (x<0) {x=x*-1; sign=-1;}
 i=x;
 x=1;
 while (i) x=x*i--;
 return sign*x; 

}
/***/
/*** Function to summarize the MCMC ***/
void MCMC_summarize (istream& infile, int burnin) {
    char temp[100]; //to store what is read from file
    double sum_LH=0; //the sum of the likelihood values (LH)
    int number_LH=0; //the number of LH values read
    int n_substitute = 0; //the number of substitute taxa
    int n_generations = 0; //the number of generations of the chain
    char generations_err = 'F'; //if the given number of generations differes from the read number warn user
    int sampel_frequency = 0; //store sample frequency
    int iteration=0; //store which generation is being read
    data_tree netdiv; //data_tree to store net diversification values
    char netdiv_err = 'F'; //flag for errors in reading net diversification rate
    data_tree ext_prop; //data_tree to stor extinction proportion values
    char ext_prop_err = 'F'; //flag for errors in reading extinction proportions
    double acceptance_sum; //the sum of the final acceptance rates for all chains
    double acceptance_n; //the number of acceptance rates read
    double acceptance_min=1; //the minimum acceptance rate for any tree
    double acceptance_max=0; //the highest acceptance rate for any tree
    double temp_acceptence; //temp to store the acceptance rate at the present iteration
    char acceptance_err = 'F'; //flag to store if something is strange with the acceptance rate values read
    double node_acceptance_sum; //same as above but for node heights
    double node_acceptance_n;
    double node_acceptance_min=1;
    double node_acceptance_max=0;
    double node_temp_acceptence;
    char node_acceptance_err = 'F';
    int n_trees=1; //the number of trees read, assumed to be at least one

    while(infile) { //read until end of file
        infile >> temp; //read next word
        if (!strcmp(temp,"(nsubstitute):")) { //if matching (nsubstitute) read the number of substitute taxa
            infile >> temp; //the next word is the number of substitute taxa
            n_substitute = atoi(temp); //read it as an integer
        }
        else if (!strcmp(temp,"(ngen):")) { //if matching (ngen) read the number of generations
            infile >> temp; //the next word is number of generations
            n_generations = atoi(temp); //read it as an integer
        }
        else if (!strcmp(temp,"(sampfreq):")) { //if matching (sampfreq) read the sample frequency
                infile >> temp; //the next word is the sample frequency
                sampel_frequency = atoi(temp); //read it as an integer
        }
        else if ((!strcmp(temp,"Acceptance_rate") && n_substitute == 0) //If no substitute taxa and reaching the acceptance rate header start parsing the MCMC
                  ||  (!strcmp(temp,"Acceptance_rate_nodes") && n_substitute > 0)) { //or if there are substitute taxa wait until the acceptance rate for node heights is reached
            int collumn=1; //we start with the first collumn
            while (infile) {
                infile>>temp; //continue pars the file
                //Null the column count when pasing the final column
                if ( n_substitute == 0 && collumn >8 ) { collumn=1; }
                else if ( n_substitute > 0 && collumn >9 ) { collumn=1; }
                //Pars the result
                if ( collumn == 1 ) { //if at the first collumn pars which tree we are at
                    if (atoi(temp)>n_trees) { //if we have reached a new tree
                        n_trees = atoi(temp); //store the number
                        acceptance_sum += temp_acceptence; //new tree time to summarize the acceptance rate for the previous tree
                        ++acceptance_n; //keep track of the number (feels safer than using n_trees)
                        if (temp_acceptence<acceptance_min) acceptance_min = temp_acceptence; //if acceptance rate is the lowest so far store it
                        if (temp_acceptence>acceptance_max) acceptance_max = temp_acceptence; //if it is the highest store it
                        if (n_substitute > 0) { //if there are substitute taxa do the same for their acceptance rate
                            node_acceptance_sum += node_temp_acceptence;
                            ++node_acceptance_n;
                            if (node_temp_acceptence<node_acceptance_min) node_acceptance_min=node_temp_acceptence;
                            if (node_temp_acceptence>node_acceptance_max) node_acceptance_max=node_temp_acceptence;
                        }
                    }
                }
                else if ( collumn == 2 ) { //if at the second collumn parse which iteration we are at
                    iteration = atoi(temp);
                    if ( iteration > n_generations || iteration < 0 ) generations_err = 'T'; //if we reach a iteration number higher than the number of generations something is wrong (or negative value)
                }
                else if ( collumn == 3 && iteration >  burnin ) { sum_LH += 1/atof(temp); ++number_LH; } //if column 3 add to the sum of LH values, and the number of LH values read
                else if ( collumn == 6 && iteration >  burnin ) { //if column 6 read net diversification rate
                    netdiv.insert_data(atof(temp)); //add value to the netdiv data_tree
                    if (atof(temp) < 0) netdiv_err = 'T'; //if negative value issue warning
                }
                else if ( collumn == 7 && iteration >  burnin ) { //if column 7 read extinction proportion
                    ext_prop.insert_data(atof(temp)); //add value to extinction proportion tree
                    if (atof(temp) < 0 || atof(temp) > 1) ext_prop_err = 'T'; //if extinction proportion is over 1 or negative issue warning
                }
                else if ( collumn == 8 && iteration >  burnin ) { //if column 8 read acceptance rate
                    temp_acceptence = atof(temp); //store it temporarily
                    if (temp_acceptence < 0 || temp_acceptence > 1) acceptance_err = 'T'; //if acceptance rate is more than 1 or negative issue warning
                }
                else if ( collumn == 9 && iteration >  burnin ) { //if column 9 (should only be reached if there is substitute taxa) read acceptance rate for node heights
                    node_temp_acceptence = atof(temp); //temporarily store the acceptance rate
                    if ( node_temp_acceptence < 0 || node_temp_acceptence > 1 ) node_acceptance_err = 'T'; //if acceptance rate is more than 1 or negative issue warning
                }
                ++collumn; //next we will do next collumn
            }
        }
    }
    acceptance_sum += temp_acceptence; //must summarize acceptance rate for last tree
    ++acceptance_n; //one more acceptance rate added
    if (temp_acceptence<acceptance_min) acceptance_min=temp_acceptence; //is it the lowest value so far?
    if (temp_acceptence>acceptance_max) acceptance_max=temp_acceptence; //is it the highest value so far?
    if (n_substitute > 0) { //is there are substitute taxa do it for node heights too
        node_acceptance_sum += node_temp_acceptence;
        ++node_acceptance_n;
        if (node_temp_acceptence<node_acceptance_min) node_acceptance_min=node_temp_acceptence;
        if (node_temp_acceptence>node_acceptance_max) node_acceptance_max=node_temp_acceptence;
    }
    //start summarizing
    cout << "Summarizing results" << endl << "********************************" << endl;
    if ( netdiv_err == 'T' ) cout << "WARNING!!! Net diversification rate had strange values!!!" << endl << endl; //if net diversification warning issued print it
    if ( ext_prop_err == 'T' ) cout << "WARNING!!! Relative extinction rate had strange values!!!" << endl << endl; //if extinction proportion warning issued print it
    if ( generations_err == 'T' ) cout << "WARNING!!! The number of generations read does not agree with what was given at start!!!" << endl << endl; //if warning issued for the generation number read print it
    if ( acceptance_err == 'T' ) cout << "WARNING!!! The acceptance rate had strange values!!!" << endl << endl; //if the acceptance rate had strange values print a warning
    if ( node_acceptance_err == 'T' ) cout << "WARNING!!! The acceptance rate for node heights had strange values!!!" << endl << endl; //same as above for node heights

    cout << "The MCMC to summarize had " << n_generations << " generations, and were sampled every " << sampel_frequency << " generations." << endl; //basic statistics on the chain
    if (n_trees == 1) cout << "The given values are summarized over all trees." << endl; //good to know, right?
    cout << "Harmonic mean of the log likelihood: " << number_LH/sum_LH << endl; //read the text
    if (n_trees == 1) cout << "The acceptance rate: " << acceptance_sum << endl; //if only one tree
    else if (n_trees > 1) { //if more than one tree
        cout << "Number of trees: " << n_trees << endl;
        cout << "Mean acceptance rate for all trees: " << acceptance_sum/acceptance_n << " (min: " << acceptance_min << ", max:" << acceptance_max << ")" << endl;
    }
    if ( n_substitute > 0 ) { //if substitute taxa give statistics for their estimation
            cout << "Number of substitute taxa: " << n_substitute << endl;
            if (n_trees == 1) cout << "Acceptance rate for node heights: " << node_acceptance_sum << endl; //if only one tree
            else if (n_trees > 1) { //if more than one tree to summarize
            cout << "Number of trees: " << n_trees << endl;
            cout << "Mean acceptance rate for node heights for all trees: " << node_acceptance_sum/node_acceptance_n << " (min: " << node_acceptance_min << ", max:" << node_acceptance_max << ")" << endl;
        }
    }
    //read the text
    cout << "Number of values read for net diversification estimation: " << netdiv.n_nodes() << ", and for relative extinction rate: " << ext_prop.n_nodes() << endl;
    cout << "Mean net diversification rate: " << netdiv.mean() << ", median: " << netdiv.median() << ", standard deviation: " << netdiv.sd() << ", 95\% credibility interval: " << netdiv.get_value(int(netdiv.n_nodes()*0.025)) << "-" << netdiv.get_value(int(netdiv.n_nodes()*0.975)) << ", min: " << netdiv.min() << ", max: " << netdiv.max() << endl;
    cout << "Mean relative extinction rate: " << ext_prop.mean() << ", median: " << ext_prop.median() << ", standard deviation: " << ext_prop.sd() << ", 95\% credibility interval: " << ext_prop.get_value(int(ext_prop.n_nodes()*0.025)) << "-" << ext_prop.get_value(int(ext_prop.n_nodes()*0.975)) << ", min: " << ext_prop.min() << ", max: " << ext_prop.max() << endl;

    cout << endl;

   //over and out
}

/*** Function to print help message ***/
void help() {

   cout << "This is SubT version 1.1 released June 2011."<< endl << endl;
   cout << "usage: SubT [file]" << endl << endl;
   cout << "SubT will read commands and data from an infile. The file must contain the number" << endl;
   cout << "   of taxa there are node height data for (e.g. ntaxa 10) and data." << endl;
   cout << "The command start will tell SubT to start reading the data and initiate the MCMC analysis." << endl;
   cout << "The data should be space-separated and there should be ntaxa -1 node heights for each set of" << endl;
   cout << "   node heights (i.e. for node heights from separate trees). Alternatively newick formated can" << endl;
   cout << "   be given as indata." << endl;
   cout << "The command end will tell the program to stop reading data and do no further MCMC." << endl;
   cout << "   (e.g. start" << endl;
   cout << "         10.88 10.72 10.50 7.36 4.86 2.02 1.38 0.26 0.000000001" << endl;
   cout << "         end   )" << endl;
   cout << endl;
   cout << "Other commands are:" << endl;
   cout << "ngen - the number of generations to run each MCMC (e.g. ngen 1000000)." << endl;
   cout << "sampfreq - the frequency with which to sample the MCMC (e.g. sampfreq 1000)." << endl;
   cout << "model - use the Yule (y) or Birth-Death (d) model (e.g. model d)." << endl;
   cout << "nsubstitute - the number of taxa in the clade for which there are no node height data (e.g." << endl;
   cout << "   nsubstitute 5)." << endl;
   cout << "windowidth - the multiplier for the proposal distribution for the net diversification" << endl;
   cout << "   rate and relative extinction rate (e.g. windowidth 0.01)." << endl;
   cout << "nodewindow - the multiplier for the proposal distribution for the node heights of the" << endl;
   cout << "   substitute taxa (e.g. nodewindow 0.5)." << endl;
   cout << "meanGamma - the mean of the net diversification rate prior distribution (e.g. meanGamma 0.1)." << endl;
   cout << "shapeGamma - the shape (as an integer) of the net diversification rate prior distribution" << endl;
   cout << "   (e.g. shapeGamma 3)." << endl;
   cout << "nodelog - if defined, the file in which to store the node heights of substitute taxa (e.g. nodelog node.log)." << endl;
   cout << endl;

}

/***********************************
*** Funktions for thr tree class ***
***********************************/
/** Function to delete tree ***
*** begining with leaf      **/
void tree::destroy_tree (node *leaf) {
    if (leaf != 0 ) {
        destroy_tree(leaf->left); //Recursively destroy left nodes
        destroy_tree(leaf->right); //Recursively destroy right nodes
        delete leaf; //Acctually destroy present node
    }
}

/** Function to get nodeheights ***
*** The function is destructive ***
*** deleting the tree in the    ***
*** process                     **/
double *tree::node_heights ( double nodeheights[], int n_nodes ) {
    double branchlength; //variable to store the sum of branchlength down to present node
    node *present = root; //store present node in tree, start with root
    node *right; //store the right node location
    int i =0; //counter

    while ( !(present->left == 0 && present->right == 0 && present->parent == 0) ) { //while not at the root with both leaf deleted
        if (present->left == 0 && present->right == 0) { //if at a tip or previously visited node
            branchlength = present->branchlength; // store the branchlength
            node *temp = present; //store present location
            present = present->parent; //move back a node
            delete temp; //delete the node just visited
            if (!(present->left == 0 && present->right == 0)) { //if we have recorded this node before
                right = present->right; //store which node is to the right
                present->left = 0; //null the left node
                present->right = 0; // and null the right node to signal that we have been here before
                present->branchlength += branchlength; //summerize the length to the tip
                if ( i > n_nodes) { cerr << "Function node_heights in tree (class) tried to write outside allocated memory; quitting without finishing reading node heights." << endl; return nodeheights; }; //make sure that we will not write outside the allocated memory
                nodeheights[i++] = branchlength; //store the node height of the node
                present = right; //move to the node to the right
            }
        }
        else { while ( present->left != 0 ) {present = present->left; } } // if not at tip or previously visited node go to the left most tip
    }
    sort(nodeheights,nodeheights+n_nodes); //sort the nodetips
    reverse(nodeheights,nodeheights+n_nodes); //reverse the order so the oldest tip is first
    return nodeheights;
}

/** Function to parse tree from newik formated file                   ***
*** Basically follow a algorithm given by Revell at:                  ***
*** http://phytools.blogspot.com/2011/02/building-tree-in-memory.html **/
void tree::tree_file_parser( istream& infile ) {
    node *present_node = root; //the tree will be parsed from the root
    node *parent_node; //keep track of wich node we came from
    char temp; //store input one character at the time
    infile >> temp; //read input one character at the time
    while (infile) {
        //ignore things in square brachets
        if ( temp == '[' ) {
            int n_right_square_brachets = 1; //keep track of number of right brachets
            while (infile && n_right_square_brachets > 0) { //while more right than left brachets have been read
                infile >> temp; //keep reading input one character at the time
                if ( temp == '[' ) ++n_right_square_brachets; //count right brachets
                if ( temp == ']' ) --n_right_square_brachets; //a left brachet neutralizes a right brachet
            }
        }
        //if left brachet create a new left node
        else if ( temp == '(' ) {
            present_node->left = new node; //create the left node
            parent_node = present_node; //store present node
            present_node = present_node->left; //move to new node
            present_node->parent = parent_node; //set previous node as parent node
            present_node->left = 0; //give left node null pointer
            present_node->right = 0; //give right node null pointer
        }
        //if , go back one node and create new right node
        else if ( temp == ',' ) {
            present_node = present_node->parent; //go back one node
            present_node->right = new node; //create new right node
            parent_node = present_node; //store present node location
            present_node = present_node->right; //go to newly created node
            present_node->parent = parent_node; //give previous node as parent node
            present_node->left = 0; //set null pointer for left node
            present_node->right = 0; //set null pointer for right node
        }
        //if right brachet move back to parent node
        else if ( temp == ')' ) {
            present_node = present_node->parent; //move to parent node
        }
        //if colon read branchlength
        else if ( temp == ':' ) {
            int i=0; //counter
            char number[100]; //branchlength can be max 99 digits
                infile >> temp; //read next character
            while ( (temp == '0' || temp == '1' || temp == '2' || temp == '3' || temp == '4' ||
                    temp == '5' || temp == '6' || temp == '7' || temp == '8' || temp == '9' ||
                    temp == '.') && i < 99 ) { //if a digit store it up to 99 digits
                number[i] = temp; //store digit
                number[++i] = '\0'; //increes counter and set next slot to end of string
                infile >> temp; //read next character
            }
            present_node->branchlength = atof(number); //when done reading branchlength in store it
            continue; //we have already read the next character so skipp that
        }
        //if we read semicolon we have reached the end of the tree
        else if ( temp == ';' ) {
            if (present_node == root ) {break;} //check so we are back at root, if so we happilly finish
            else { //if we are not back at the root something is wrong
                cerr << "Right and left brackets not balanced in tree file, removing tree." << endl; //print error measeage
                destroy_tree(root); //destroy what we have built
                break; //finish this shit
            }
        }
        infile >> temp; //read next character for next loop
    }
}

