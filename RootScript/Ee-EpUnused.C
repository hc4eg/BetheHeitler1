  //vector<Monitor>* moni0 = new vector<Monitor>;
  //vector<Monitor>* moni1 = new vector<Monitor>;
  
  //  BH_Event* event = 0;

  //TBranch* bmoni0 = 0;
  //TBranch* bmoni1 = 0;

    cerr << "Event_Branch is set." << endl;
  //  bmoni0 = tree -> GetBranch("monitor0");
  //  bmoni1 = tree -> GetBranch("monitor1");
  //cerr << "monitors branch is set, with address: " << &moni0 << endl;

    //    bmoni0 -> SetAddress(&moni0);
    //    bmoni1 -> SetAddress(&moni1);

      cerr << "mE0 = " << mE[0] << endl;
      cerr << "mE1 = " << mE[1] << endl;
      cerr << "moni0E = " << moni0->energy_m << endl;
      cerr << "moni1E = " << moni1->energy_m << endl;

    //      delete moni0; moni0 = new Monitor;
    //      delete moni1; moni1 = new Monitor;

    //    moni0->ClearMonitor();
    //    moni1->ClearMonitor();

    //          bmoni0 -> GetEvent(i);
    //          bmoni1 -> GetEvent(j);
      //cerr << "Reading the " << i << "th entry." << endl;
      //cerr << "Size of monitor at this event is " << event->monitors.size() << endl;

      //      mcharge[0] = 0; mcharge[1] = 0;

      /*
	  if(moni0->enumber == moni1->enumber) 
	    {
	      mE[0] = moni0->energy_m;
	      mE[1] = moni1->energy_m;
	    }
	  else if(moni0->enumber > moni1->enumber) 
	    do{
	        j++;
		bmoni1 -> GetEntry(j);
	    }while(moni0->enumber != moni1->enumber);
	  else 
	    do{
	        i++;
		bmoni0 -> GetEntry(i);
	    }while(moni0->enumber != moni1->enumber);
      */



      //     cerr << "Processing " << i+1 << "th event." << endl;
      //     cerr << "Ep = " << mE[0] << endl;
      //     cerr << "Ee = " << mE[1] << endl;

      //cerr << << "Charge in monitor is " << monit0->at(0).charge << endl;
      
      /*  Int_t j = 0;
      for (vector<Monitor>::iterator it = moni0->begin() ; it != moni0->end(); ++it)
    {
      mcharge[j] = (*it)->charge;
      mE[j] = (*it)->energy_m;
      j++;
    }
      
      */

      /*
      for (unsigned j = 0 ; j < monitor->size() ; j++)
    {
      mcharge[j] = monitor->at(j).charge;
      mE[j] = monitor->at(j).energy_m;
      cerr << "Size of monitor at this event is " << monitor->size() << endl;
      cerr << "Charge in monitor is " << mcharge[j] << endl;
      cerr << "Energy is " << mE[j] << endl;
      j++;
    }
      */
      
      //if((mcharge[0] != 0) && (mcharge[1] != 0))
