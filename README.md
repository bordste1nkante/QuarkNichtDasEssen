# QuarkNichtDasEssen
# hier könnten ihre Kommentare stehen


Grundlegende Struktur (Mein Vorschlag):

# This sketch is with metropolis algorithm in mind
It is likely resourceful, if we do not generate an X for each update, but an array of matrices from which we uniformly pick on. After n sweeps we update parts or the whole list.

void main()

params = read_out_yamlOrTxt(); //pass a file that contains all important information like lattice size and where to save data

lattice = array of current lattice state //have to figure out how precisly it should be done, maybe array of array pointers, which point to matrices

//depending on params, either cold or hot start
cold_start(); or hot_start();

for n in numberOfSweeps{
    for i in latticeSites //likely multiple for loops, one per dimension
    X = random_SU(3)_generator(); 
    U' = XU; // U is previous matrix at lattice site (technically links between lattice sites/physically)
    bool accept = accept_reject(i,U',&lattice); //return true (accepted) or false (rejected)
    if bool == true{
        lattice[i] = U'
    }
}
//might add an additional for loop to account for sweeps or increase numberOfGeneratedConfigurations and use modulo x to extract only every x-th
for n in numberOfGeneratedConfigurations{
    for i in latticeSites //likely multiple for loops, one per dimension
    X = random_SU(3)_generator(); 
    U' = XU; // U is previous matrix at lattice site (technically links between lattice sites/physically)
    bool accept = accept_reject(i,U', &lattice); //return true (accepted) or false (rejected)
    if bool == true{
        lattice[i] = U';
        store_data(U');

    }
    else{
        store_data(U)
    }
return;
}


//sets all matrices to identity
void cold_start(&lattice){
    for i in latticeSites{
        lattice[i]=identityMatrix;
    }
    return;

}

//sets all matrices to randomSU(3)
void hot_start(&lattice){
    for i in latticeSites{
        lattice[i]=random_SU(3)_generator();
    }
    return;

}
//generates a randomSU(3) and its inverse
someVectorShit random_SU(3)_generator(){
    Y = some algorithm;
    Y^-1 = some algorithm;
    X =select uniformly from Y and Y^-1;
    return X;

}

//accepts the new configuration, if action decreases or with some chance
bool accept_reject(i, U', const &lattice){
    U = lattice[i];
    SOfU = calculateAction(U, lattice,i);
    SOfU' = calculateAction(U', lattice,i);
    
    accept = true if condition met else false;
    return accept;

}


//calculates action, see Chromodynamics on the lattice for lattice action in terms of U's
double calculateAction(U,lattice,i){
    S= calculation
    return S
}








