cell_array::cell_array(int Nxx, int Ntypess, cell &c) : next_diffusion_events(vector<diffusionstore2>(Nxx*Nxx)),reactors(matrix<cell>(Nxx, Nxx)), si(vector1<self_interaction *>(Ntypess)), cross_si(matrix<self_interaction *>(Ntypess,Ntypess)), cnni(matrix<self_interaction *>(Ntypess,Ntypess))
{
    Nx=Nxx;
    Ntypes = Ntypess;
    current_time = 0.;
    for (int i = 0; i < Nxx; i++)
    {
        for (int j = 0; j < Nxx; j++)
        {
            reactors(i, j) = c;
        }
    }

    for (int i = 0; i < Ntypes; i++)
    {
        no_self_attraction *pot1 = new no_self_attraction;
        si[i] = pot1->clone();
        delete pot1;
        for (int j = 0; j < Ntypes; j++)
        {
            no_surface_interaction *pot2 = new no_surface_interaction;
            cross_si(i, j) = pot2->clone();
            cnni(i, j) = pot2->clone();
            delete pot2;
        }
    }
};

void cell_array::initialize_next_diffusion_events() {
    
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Nx; j++)
        {
            diffusionstore2 sd2;
            sd2.posi = i;
            sd2.posj = j;
            //we now store the position at which it occurs
            reactors(i, j).sorted_times[0];
            
            sd2.a = reactors(i,j).sorted_times[0];
            
            next_diffusion_events[i * Nx + j] = sd2;
            //next_diffusion_events[i*Nx+j].a = reactors(i,j).sorted_times[0]; //add the minimal time for all the reactors
        }
    }
    std::sort(next_diffusion_events.begin(),next_diffusion_events.end() ); //sort all these times
}


int periodic_correct(int i, int N) {
    if (i<0) return i+N;
    else if(i > N-1) return i-N;
    else {return i;}

}

int checkEquality(int i1, int i2, int j1, int j2)
{
    return (i1 == i2 && j1 == j2) ? 1 : 0;
}

double cell_array::getenergy(int i1, int j1, int i2, int j2, int tm) {

    //i1,j1 are the place the particles was, (i2,j2) are where the particle will move to
    //we subtract one index of type tm from i1,j1 and add one index of type tm to (i2,j2)

    int indexi = i1;
    int indexj = j1;

    int index1nn1i =  periodic_correct(i1 - 1,Nx);
    int index1nn1j = j1;

    int index1nn2i = periodic_correct(i1 + 1,Nx);
    int index1nn2j = j1;

    int index1nn3i = i1;
    int index1nn3j = periodic_correct(j1 - 1,Nx);

    int index1nn4i = i1;
    int index1nn4j = periodic_correct(j1 + 1,Nx);

    int index2nn1i = periodic_correct(i2 - 1, Nx);
    int index2nn1j = j2;

    int index2nn2i = periodic_correct(i2 + 1, Nx);
    int index2nn2j = j2;

    int index2nn3i = i2;
    int index2nn3j = periodic_correct(j2 - 1, Nx);

    int index2nn4i = i2;
    int index2nn4j = periodic_correct(j2 + 1, Nx);

    double befenergy = 0.0;
    double aftenergy = 0.0;
    //we only need to consider the type that moves (all other volume energies remain the same)

    befenergy += (*si[tm])(reactors(i1, j1).typenumber[tm]); //volume interaction, self
    befenergy += (*si[tm])(reactors(i2, j2).typenumber[tm]); // volume interaction, self

    aftenergy += (*si[tm])(reactors(i1, j1).typenumber[tm] - 1); // volume interaction, self
    aftenergy += (*si[tm])(reactors(i2, j2).typenumber[tm] + 1); // volume interaction, self


    for(int i = 0  ; i < Ntypes ; i++) {
        if(i!=tm) {
            befenergy += (*cross_si(i, tm))(reactors(i1, j1).typenumber[i], reactors(i1, j1).typenumber[tm]);      // volume interaction, other
            befenergy += (*cross_si(i, tm))(reactors(i2, j2).typenumber[i], reactors(i2, j2).typenumber[tm]);      // volume interaction, other

            aftenergy += (*cross_si(i, tm))(reactors(i1, j1).typenumber[i], reactors(i1, j1).typenumber[tm] - 1);           // volume interaction, other
            aftenergy += (*cross_si(i, tm))(reactors(i2, j2).typenumber[i], reactors(i2, j2).typenumber[tm] + 1);           // volume interaction, other
        }
    }

    

    for(int i = 0 ; i < Ntypes ; i++) {
        
            //surface interactions

            // befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn1i, index1nn1j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn2i, index1nn2j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn3i, index1nn3j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn4i, index1nn4j).typenumber[tm]);

            // befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn1i, index2nn1j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn2i, index2nn2j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn3i, index2nn3j).typenumber[tm]);
            // befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn4i, index2nn4j).typenumber[tm]);
            befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm], reactors(index1nn1i, index1nn1j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm], reactors(index1nn2i, index1nn2j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm], reactors(index1nn3i, index1nn3j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm], reactors(index1nn4i, index1nn4j).typenumber[i]);

            befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm], reactors(index2nn1i, index2nn1j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm], reactors(index2nn2i, index2nn2j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm], reactors(index2nn3i, index2nn3j).typenumber[i]);
            befenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm], reactors(index2nn4i, index2nn4j).typenumber[i]);

            // aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn1i, index1nn1j).typenumber[tm] - 1);
            // aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn2i, index1nn2j).typenumber[tm] - 1);
            // aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn3i, index1nn3j).typenumber[tm] - 1);
            // aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[i], reactors(index1nn4i, index1nn4j).typenumber[tm] - 1);

            // aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn1i, index2nn1j).typenumber[tm] + 1);
            // aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn2i, index2nn2j).typenumber[tm] + 1);
            // aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn3i, index2nn3j).typenumber[tm] + 1);
            // aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[i], reactors(index2nn4i, index2nn4j).typenumber[tm] + 1);
            if(i == tm) { //accounts for the change in state by moving
                
                aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, checkEquality(i2, index1nn1i, j2, index1nn1j) + reactors(index1nn1i, index1nn1j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, checkEquality(i2, index1nn2i, j2, index1nn2j) + reactors(index1nn2i, index1nn2j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, checkEquality(i2, index1nn3i, j2, index1nn3j) + reactors(index1nn3i, index1nn3j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, checkEquality(i2, index1nn4i, j2, index1nn4j) + reactors(index1nn4i, index1nn4j).typenumber[i]);

                aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, -checkEquality(i1, index2nn1i, j1, index2nn1j) + reactors(index2nn1i, index2nn1j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, -checkEquality(i1, index2nn2i, j1, index2nn2j) + reactors(index2nn2i, index2nn2j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, -checkEquality(i1, index2nn3i, j1, index2nn3j) + reactors(index2nn3i, index2nn3j).typenumber[i]);
                aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, -checkEquality(i1, index2nn4i, j1, index2nn4j) + reactors(index2nn4i, index2nn4j).typenumber[i]);
            }
            
            else{
            aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, reactors(index1nn1i, index1nn1j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, reactors(index1nn2i, index1nn2j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, reactors(index1nn3i, index1nn3j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i1, j1).typenumber[tm] - 1, reactors(index1nn4i, index1nn4j).typenumber[i]);

            aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, reactors(index2nn1i, index2nn1j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, reactors(index2nn2i, index2nn2j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, reactors(index2nn3i, index2nn3j).typenumber[i]);
            aftenergy += (*cnni(i, tm))(reactors(i2, j2).typenumber[tm] + 1, reactors(index2nn4i, index2nn4j).typenumber[i]);
            }
            //we are missing one contribution to the energy because one of the nearest neighbors has changed as well
    }

    return aftenergy-befenergy;


}



void cell_array::find_next_reaction_event()
{ //can we regenerate the entire thing from new, and thus not have to loop through the entire code again?
    // int i1 = 0;
    // int j1 = 0;

    // double rmin = reactors(i1, j1).reaction_clock;

    // for (int i = 0; i < Nx; i++)
    // {
    //     for (int j = 0; j < Nx; j++)
    //     {
            
 
    //         if (reactors(i, j).reaction_clock < rmin)
    //         {
    //             rmin =  reactors(i,j).reaction_clock;
    //             i1 = i;
    //             j1 = j;
    //         }
    //     }
    // }

    // pp res;
    // res.i = i1;
    // res.j = j1;
    // res.time =  rmin;

    // next_reaction_event = res;

}

void cell_array::find_next_diffusion_event()
{
    // int i1 = 0;
    // int j1 = 0;

    // double dmin = reactors(i1, j1).diffusion_clock;
    // for (int i = 0; i < Nx; i++)
    // {
    //     for (int j = 0; j < Nx; j++)
    //     {

    //         if (reactors(i, j).diffusion_clock < dmin)
    //         {
    //             dmin = reactors(i, j).diffusion_clock;
    //             i1 = i;
    //             j1 = j;
    //         }
    //     }
    // }

    // next_events[0].
    // pp res;
    // res.i = i1;
    // res.j = j1;
    // res.time = dmin;

    // next_diffusion_event = res;
}

void cell_array::diffusionstay(int indexi,int indexj,int index_type_to_move) {
 //the particle stays
    
    //we erase the things associated with the old particle
    reactors(indexi, indexj).sorted_times.erase(reactors(indexi, indexj).sorted_times.begin());
    next_diffusion_events.erase(next_diffusion_events.begin()); // we erase it from this too;

    //we generate a new time
    double time1 = reactors(indexi, indexj).sd[index_type_to_move].generate_alarm();
    diffusionstore sd;
    sd.nt = index_type_to_move;
    sd.time = current_time + time1; 

    //we add the new time to the list
    insert_new(reactors(indexi, indexj).sorted_times, sd);

    //we insert the new back to the list

    diffusionstore2 sd3;
    sd3.posi = indexi;
    sd3.posj = indexj;
    sd3.a = reactors(indexi,indexj).sorted_times[0];

    
    // cout << "added over empty stay " << sd3.a.nt << endl;


    auto insertion_point = std::lower_bound(next_diffusion_events.begin(), next_diffusion_events.end(), sd3);
    next_diffusion_events.insert(insertion_point, sd3);
    
}

void cell_array::diffusionmove(int indexi,int indexj, int indexmovei, int indexmovej, int index_type_to_move) {

     //WE ARE NOT ADDING BACK THE DIFFUSION WHEN IT STAYS IN PLACE FOR SOME REASON


    //do a diffusive move from i1,j1 to i2,j2

    //first thing to do is remove the elements from the list now the move has been done
    // cout << reactors(indexi, indexj).sorted_times.size() << " " << reactors(indexmovei, indexmovej).sorted_times.size() << endl;

    // cout << reactors(indexi,indexj).sorted_times.size() << " " << reactors(indexmovei,indexmovej).sorted_times.size()  << endl;

    reactors(indexi, indexj).sorted_times.erase(reactors(indexi, indexj).sorted_times.begin());
    next_diffusion_events.erase(next_diffusion_events.begin()); // we erase it from this too;


    //the numbers in each of the boxes have to be updated
    reactors(indexi, indexj).typenumber[index_type_to_move]--;
    reactors(indexmovei, indexmovej).typenumber[index_type_to_move]++;

    //next we set find out what the next time for the motion of the diffused particle is
    double time1 = reactors(indexmovei, indexmovej).sd[index_type_to_move].generate_alarm();
    diffusionstore sd;
    sd.nt = index_type_to_move;
    sd.time = current_time + time1;

    //we now need to update both lists, the thing to check is if the time for indexmovei,indexmovej 
    // is less than the prior time
    if(reactors(indexmovei,indexmovej).sorted_times.size() > 0 ) { //if the moving list exists
    if(sd.time < reactors(indexmovei,indexmovej).sorted_times[0].time) {
        //we need to erase the previous jump at indexmovei,indexmovej if the new time is less
            int iter;
            bool cond1=false;
            bool cond2=false;
            for(int i = 0 ; i < next_diffusion_events.size() ; i++) {
                cond1 = next_diffusion_events[i].posi == indexmovei;
                cond2 = next_diffusion_events[i].posj == indexmovej;
                if(cond1 && cond2) {
                    iter = i;
                    break;
                }
            }
            //only if there is a match do we do thei following
            if(cond1 && cond2) {
            next_diffusion_events.erase(next_diffusion_events.begin() + iter); //erase the old time
            //add the new time
            diffusionstore2 sd2;
            sd2.posi = indexmovei;
            sd2.posj = indexmovej;
            sd2.a = sd;

            // cout << "added over previous jumped " << sd2.a.nt << endl;

            
            auto insertion_point = std::lower_bound(next_diffusion_events.begin(), next_diffusion_events.end(), sd2);
            // cout << insertion_point - next_diffusion_events.begin() << endl;
            // pausel();
            next_diffusion_events.insert(insertion_point, sd2);

            }
        }      
    }
    else{ //if it doesn't add anyway

        diffusionstore2 sd2;
        sd2.posi = indexmovei;
        sd2.posj = indexmovej;
        sd2.a = sd;

        // cout << "added over empty jumped" << sd2.a.nt << endl;

        auto insertion_point = std::lower_bound(next_diffusion_events.begin(), next_diffusion_events.end(), sd2);

        next_diffusion_events.insert(insertion_point, sd2);

    }


    //add the new time to the local information in the reactor
    insert_new(reactors(indexmovei, indexmovej).sorted_times, sd);

    
    // cout << reactors(indexi,indexj).sorted_times.size() << " " << reactors(indexmovei,indexmovej).sorted_times.size()  << endl;
    

    // now we need to add the new diffusion event corresponding to posi,posj (where the particle moved from)


        if(reactors(indexi,indexj).sorted_times.size() > 0) { //only add back if it's non empty
        //cout << "added back" << endl;
        diffusionstore2 sd3;
        sd3.posi = indexi;
        sd3.posj = indexj;
        sd3.a = reactors(indexi,indexj).sorted_times[0];

        
        // cout << "added over empty stay " << sd3.a.nt << endl;


        auto insertion_point = std::lower_bound(next_diffusion_events.begin(), next_diffusion_events.end(), sd3);
        next_diffusion_events.insert(insertion_point, sd3);

        }
    



    // if(numb != next_diffusion_events.size() ) {
    //     pausel();
    // }
    //cout << next_diffusion_events.size() << endl;
    //this->collision_check();
    //pausel();

    //what happens if there are no particles there?


}

void cell_array::collision_check() {
    int si = next_diffusion_events.size();
    ofstream myfile;
    myfile.open("check.csv");
    for(int i = 0 ; i < si ; i++) {
        myfile << next_diffusion_events[i].posi << "," << next_diffusion_events[i].posj << "," << next_diffusion_events[i].a << endl;
    }
    myfile.close();
    
    // for(int i = 0 ; i <  si-1 ; i ++) {
        
    //         int p1 = next_diffusion_events[i].posi;
    //         int q1 = next_diffusion_events[i].posj;

    //         int p2 = next_diffusion_events[i+1].posi;
    //         int q2 = next_diffusion_events[i+1].posj;

    //         if(p1==p2 && q1 == q2) {
    //             cout << "collision detected!" << endl;
    //             pausel();
    //         }
    //         else {
    //             //cout << "no collisions" << endl;
    //         }
        
    // }
}

void cell_array::do_a_diffusion_event(vector1<int> &counter) {

   int indexi = next_diffusion_events[0].posi;
   int indexj = next_diffusion_events[0].posj;
   double dt = next_diffusion_events[0].a.time - current_time;
   current_time = next_diffusion_events[0].a.time;

   int mov1 = -1 + 2 * (rand() % 2);
   int mov2 = -1 + 2 * (rand() % 2);

   int indexmovei = periodic_correct(indexi + mov1, Nx);
   int indexmovej = periodic_correct(indexj + mov2, Nx);

   int index_type_to_move = next_diffusion_events[0].a.nt;
   //cout << index_type_to_move << endl;

   double delE = this->getenergy(indexi, indexj, indexmovei, indexmovej, index_type_to_move);

   double r = (double)rand() / (double)RAND_MAX;

   if (delE <= 0)
   {
    counter[0]++;
       // do a move
       this->diffusionmove(indexi, indexj, indexmovei, indexmovej, index_type_to_move);

       // SHOULD THESE TIMES BE COMPLETELY RESET?
       // reactors(indexi, indexj).generate_reaction_time(current_time);
       // reactors(indexmovei, indexmovej).generate_reaction_time(current_time);

       // reactors(indexi, indexj).generate_diffusion_time(); // get the  next minima in diffusion times for this cell
       // reactors(indexmovei, indexmovej).generate_diffusion_time(); // get the  next minima in diffusion times for this cell

       // for (int i = 0; i < Nx; i++)
       // {
       //     for (int j = 0; j < Nx; j++)
       //     {
       //         if ((i == indexi && j == indexj) || (i == indexmovei && j == indexmovej))
       //         {
       //         }
       //         else
       //         {
       //             reactors(i, j).reaction_clock += dt; //this is equivalent to redrawing all the times
       //         }
       //     }
       // }

    }
    else if(exp(-delE) > r )
    {
        // double time1 = reactors(indexi, indexj).sd[index_type_to_move].generate_alarm();
        // reactors(indexi, indexj).diffusion_times[index_type_to_move].erase(reactors(indexi, indexj).diffusion_times[index_type_to_move].begin() + ind);
        // reactors(indexmovei, in

        counter[1]++;
        this->diffusionmove(indexi, indexj, indexmovei, indexmovej, index_type_to_move);

        //SHOULD THESE TIMES BE COMPLETELY RESET?

        // //The times should all be reset, think of a random collision model (see notebook)
        // reactors(indexi, indexj).generate_time(current_time);
        // reactors(indexmovei, indexmovej).generate_time(current_time);

        // for(int i = 0 ; i < Nx ; i++) {
        //     for(int j = 0 ; j < Nx ; j++) {
        //         if ((i == indexi && j == indexj) || (i == indexmovei && j == indexmovej))
        //         {
        //         }
        //         else{
        //         reactors(i, j).reaction_clock += dt;
        //         }
        //     }
        // }
        
        // reactors(indexi, indexj).generate_diffusion_time();         // get the  next minima in diffusion times for this cell
        // reactors(indexmovei, indexmovej).generate_diffusion_time(); // get the  next minima in diffusion times for this cell
        }
    else{

            counter[2]++;
        this->diffusionstay(indexi, indexj, index_type_to_move);
        //when these numbers are the same, it should just rewrite the part corresponding to the (indexi,indexj)

        //we now face a key assumption about what happens to the move probability after a move was attempted, do we reset it completely?
        //we here say that the number is redrawn again
        // double time1 = reactors(indexi, indexj).sd[index_type_to_move].generate_alarm();
        // reactors(indexi, indexj).diffusion_times[index_type_to_move][ind] = time1 + current_time;

        // reactors(indexi, indexj).generate_diffusion_time();

        //state did not change, do not redraw probs
    }

    //     cout << next_diffusion_events.size() << endl;
    // int numb=0;
    // int total =0;
    // for(int i = 0  ; i < Nx ; i++) {
    //   for(int j = 0 ; j < Nx ; j++) {
    //     total += reactors(i,j).sorted_times.size();
    //     if(reactors(i,j).sorted_times.size() > 0) {
    //         numb++;
    //     }
    //   }  
    // }
    // cout << numb << " " << total << endl;
    // cout << endl;


    //this->find_next_reaction_event();  // get the next reaction event
   // this->find_next_diffusion_event(); // get the next diffusion event

}

void cell_array::reset_all_times() {
    double ct = current_time;
    for(int i = 0 ; i < Nx ; i++) { //reset all events
        for(int j = 0  ; j < Nx ; j++) {
            int n =  reactors(i,j).sorted_times.size();
            for(int k = 0 ; k < n ; k++)
            reactors(i,j).sorted_times[k].time -= ct;
        }
    }

    for (int i = 0; i < next_diffusion_events.size() ; i++)
    { // reset all events
        next_diffusion_events[i].a.time -= ct;
    }


};

void cell_array::iterate() {

    //if(next_diffusion_event.time < next_reaction_event.time ) {
        // do a diffusion event
    if(true)   {
        //do_a_diffusion_event();
       
    }
    else{ //do I do the time correctly?
        
        // do a reaction event
        // int indexi = next_reaction_event.i;
        // int indexj = next_reaction_event.j;
        // current_time += next_reaction_event.time; //this time is measured relative to the current time, whereas d is not

        // reactors(indexi, indexj).do_reaction(current_time); // do the reaction



        // reactors(indexi,indexj).generate_time(current_time); 
        
        // //after a chemical reaction has occured, again we redraw all the times
        // for (int i = 0; i < Nx; i++)
        // {
        //     for (int j = 0; j < Nx; j++)
        //     {
        //         if ((i == indexi && j == indexj))
        //         {
        //         }
        //         else
        //         {
        //             reactors(i, j).reaction_clock += next_reaction_event.time;
        //         }
        //     }
        // }

        // reactors(indexi, indexj).generate_diffusion_time(); // get the  next minima in diffusion times for this cell

        // this->find_next_reaction_event(); // get the next reaction event
        // this->find_next_diffusion_event(); //get the next diffusion event

        
        }
        
}