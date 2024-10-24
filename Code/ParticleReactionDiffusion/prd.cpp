cell::cell() : typenumber(vector1<int>(2)), sorted_times(vector<diffusionstore>())
{
    Ntypes = 2;
    Nreactions = 2;
    sd = new Subdiffusion[Ntypes];
    cr = new chemical_reaction[Nreactions];
}

cell::cell(int Nt, int Nr) : typenumber(vector1<int>(Nt)), sorted_times(vector<diffusionstore>())
{
    Ntypes = Nt;
    Nreactions = Nr;
    sd = new Subdiffusion  [Ntypes];
    cr = new chemical_reaction [Nreactions];
}

cell::cell(const cell &a) : typenumber(a.typenumber), sorted_times(a.sorted_times)
{
    Ntypes = a.Ntypes;
    Nreactions = a.Nreactions;
    sd = new Subdiffusion [Ntypes];
    cr = new chemical_reaction [Nreactions];
    for(int i = 0 ; i < Ntypes ; i++) {
        sd[i] = a.sd[i];
    }
    for (int i = 0; i < Nreactions; i++)
    {
        cr[i] = a.cr[i];
    }
}

cell &cell::operator=(const cell &a) {
    typenumber =  a.typenumber;
    Ntypes = a.Ntypes;
    Nreactions = a.Nreactions;
    sorted_times = a.sorted_times;
    delete sd;
    delete cr;

    sd = new Subdiffusion [Ntypes];
    cr = new chemical_reaction [Nreactions];

    for (int i = 0; i < Ntypes; i++)
    {
        sd[i] = a.sd[i];
    }
    for (int i = 0; i < Nreactions; i++)
    {
        cr[i] = a.cr[i];
    }
    return *this;
}

void cell::initialize(vector1<int> &numb) {
    
    typenumber  = numb;
    if(numb.getsize() != Ntypes) error("didn't initialize all values");
    for(int i  = 0 ; i < Ntypes ; i++) {
        for(int j  = 0 ; j < typenumber[i] ; j++) {
            double time = sd[i].generate_initial();

            diffusionstore a;
            a.nt = i;
            a.time = time;
            sorted_times.push_back(a);

        }
    }

    std::sort(sorted_times.begin(), sorted_times.end() );

}

void cell::set_chemical_reaction(chemical_reaction a, int i) {
    if(i < 0 || i >= Ntypes) 
        error("not setting the right reaction");
    cr[i] = a;
}

void cell::set_subdiffusion(Subdiffusion a, int i)
{
    if (i < 0 || i >= Ntypes)
        error("not setting the right reaction");
    sd[i] = a;
}

void cell::generate_reaction_time(double ct) {
    vector<double> reacs;
    double running_total = 0.0;
    
    for(int i = 0 ; i < Nreactions ; i++) {
        
        int s1 = typenumber[cr[i].i1];
   
        double totalrate;
        if(cr[i].orderin==2) {
            int s2 = typenumber[cr[i].i2];
            totalrate= cr[i].rate*s1*s2;
        }
        else{
            totalrate = cr[i].rate*s1;
        }

        running_total += totalrate;
        reacs.push_back(running_total);
    }
    double r = (double)rand() / (double)RAND_MAX;
    if(running_total > 1.E-12) {
    double tr = log(1. / (1. - r)) /running_total;
    reaction_clock = ct + tr;
    double r2 = running_total*(double)rand() / (double)RAND_MAX;
    for(int i = 0  ; i < reacs.size() ; i++ ) {
        if(r2 < reacs[i]) {
            which_reaction =  i;
            break;
        }
    }
    }
    else{
        reaction_clock = 1.E10;
        which_reaction=0;
    }

}

void cell::find_random_single(int type1, int &i)
{
    bool match = false;
    int n = sorted_times.size();

    while (!match)
    {
        int cho = rand() % n;
        if (sorted_times[cho].nt == type1)
        {
            i = cho;
            match = true;
        }
    }


}

void cell::find_random_pair(int type1, int type2, int &i, int &j) {
    bool match = false;
    int n = sorted_times.size();
    
    while (!match)
    {
        int cho = rand() % n;
        if (sorted_times[cho].nt == type1)
        {
            i = cho;
            match = true;
        }
    }

    match = false;
    while (!match)
    {
        int cho = rand() % n;
        if (sorted_times[cho].nt == type2)
        {
            j = cho;
            if(j != i) match = true;
        }
    }
}

void cell::do_reaction(double ct) {
    int i1 = cr[which_reaction].i1;

    if (cr[which_reaction].orderin == 2)
    {
        int i2 = cr[which_reaction].i2;
        
        int o1 = cr[which_reaction].o1;

        if (cr[which_reaction].orderout == 2)
        {
            //2 to 2
            int o2 = cr[which_reaction].o2;
            // int iter1 = std::rand() % diffusion_times[i1].size();
            // int iter2 = std::rand() % diffusion_times[i2].size();
            
            //perhaps it's faster to find a random element in the sorted list that matches?
            int iter1,iter2;
            find_random_pair(i1,i2,iter1,iter2);
            typenumber[i1]--; // chemical reaction occured, remove 1
            typenumber[i2]--; //chemical reaction occured, remove 1
            typenumber[o1]++;
            typenumber[o2]++;

            if(iter1 < iter2) std::swap(iter1,iter2);

            sorted_times.erase(sorted_times.begin() + iter1); //we want to remove the largest iterator first

            sorted_times.erase(sorted_times.begin() + iter2); //then the smallest iterator

            //what do we do to the times that we delete?
            
            double time1 = sd[o1].generate_alarm();
            double time2 = sd[o2].generate_alarm();

            diffusionstore st1;
            st1.nt=o1;
            st1.time =  time1;



            diffusionstore st2;
            st2.nt = o2;
            st2.time = time2;


            insert_new(sorted_times, st1);
            insert_new(sorted_times, st2);

            
        }
        else{
            //2 to 1
            // int iter1 = std::rand() % diffusion_times[i1].size();
            // int iter2 = std::rand() % diffusion_times[i2].size();

            // diffusion_times[i1].erase(diffusion_times[i1].begin() + iter1);
            // diffusion_times[i2].erase(diffusion_times[i2].begin() + iter2);
            int iter1, iter2;
            find_random_pair(i1, i2, iter1, iter2);
            typenumber[i1]--; // chemical reaction occured, remove 1
            typenumber[i2]--; // chemical reaction occured, remove 1
            typenumber[o1]++;

            if (iter1 < iter2)
                std::swap(iter1, iter2);
            sorted_times.erase(sorted_times.begin() + iter1);
            sorted_times.erase(sorted_times.begin() + iter2);

            double time1 = sd[o1].generate_alarm();

            diffusionstore st1;
            st1.nt = o1;
            st1.time = time1;


            insert_new(sorted_times, st1);

        }
    }
    else{
        int o1 = cr[which_reaction].o1;

        if (cr[which_reaction].orderout == 2)
        {
            //1 to 2
            int o2 = cr[which_reaction].o2;
            // int iter1 = std::rand() % diffusion_times[i1].size();

            // diffusion_times[i1].erase(diffusion_times[i1].begin() + iter1);
            int iter1;
            find_random_single(i1, iter1);
            typenumber[i1]--;
            typenumber[o1]++;
            typenumber[o2]++;

            sorted_times.erase(sorted_times.begin() + iter1);

            double time1 = sd[o1].generate_alarm();
            double time2 = sd[o2].generate_alarm();

            diffusionstore st1;
            st1.nt = o1;
            st1.time = time1;


            diffusionstore st2;
            st2.nt = o2;
            st2.time = time2;


            insert_new(sorted_times, st1);
            insert_new(sorted_times, st2);

        }
        else
        {
            //1 to 1
            // int iter1 = std::rand() % diffusion_times[i1].size();

            // diffusion_times[i1].erase(diffusion_times[i1].begin() + iter1);
            int iter1;
            find_random_single(i1, iter1);
            typenumber[i1]--;
            typenumber[o1]++;

            sorted_times.erase(sorted_times.begin() + iter1);

            double time1 = sd[o1].generate_alarm();

            diffusionstore st1;
            st1.nt = o1;
            st1.time = time1;


            insert_new(sorted_times, st1);
        }
    }
    this->generate_reaction_time(ct);
    //generate new distribution after the reaction is done
}

void cell::generate_diffusion_time() { //this finds the minimal diffusion time inside the reactor
    // int index1;
    // int index2;
    // double dmin;
    // int null_vecs = 0;
    // for (int i = 0; i < Ntypes; i++)
    // {
    //     if( diffusion_times[i].size() > 0 ) {
    //     dmin = diffusion_times[i][0];

    //     index1 = i;
    //     index2 = 0;
    //     }
    //     else{
    //         null_vecs++;
    //     }
    // }
    // if(null_vecs == Ntypes) {
    //     diffusion_clock = 1.E10;
    //     which_diffusion_type = 0;
    //     which_diffusion = 0;
    // }
    // else{
    // //what to do if the initial vector doesn't necessarily exist?

    // for (int i = 0 ; i < diffusion_times.getsize(); i++)
    // {
    //     for (int j = 0; j < diffusion_times[i].size(); j++)
    //     {
    //         if (diffusion_times[i][j] < dmin)
    //         {
    //             dmin = diffusion_times[i][j];
    //             index1 = i;
    //             index2 = j;
    //         }
    //     }
    // }
    // diffusion_clock = dmin;
    // which_diffusion_type = index1;
    // which_diffusion = index2;
    // }
}


