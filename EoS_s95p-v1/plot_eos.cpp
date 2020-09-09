  




while (!mFile.eof()) {
    mFile.getline(buff,200);
    if (!(*buff) || (*buff == '#'))
      continue;
    iss = new istringstream(buff);   
    (*iss) >> name >> mass >> gamma >> spin >> I >> I3 >> Nq >> Ns >> Naq >> Nas >> Nc >> Nac >> MC;
    number++;
    PRINT_DEBUG_2('\t'<<number<<" "<<name<<" "<<mass<<" "<<gamma<<" "<<spin<<" "<<I<<" "<<I3<<" "<<Nq<<" "<<Naq<<" "<<Ns<<" "<<Nas<<" "<<Nc<<" "<<Nac<<" "<<MC);
    tPartBuf = new ParticleType();
    tPartBuf->SetNumber(number);
    tPartBuf->SetName(name);
    tPartBuf->SetMass(mass);
    tPartBuf->SetGamma(gamma);
    tPartBuf->SetSpin(spin);
    tPartBuf->SetBarionN(static_cast<int> ((Nq + Ns + Nc)/3. - (Naq + Nas + Nac)/3.) );
    tPartBuf->SetI(I);
    tPartBuf->SetI3(I3);
    tPartBuf->SetStrangeN(static_cast<int> (Nas - Ns));
    tPartBuf->SetCharmN(static_cast<int> (Nc - Nac));
    tPartBuf->SetNumberQ(static_cast<int> (Nq));
    tPartBuf->SetNumberQ(static_cast<int> (Naq));
    tPartBuf->SetNumberQ(static_cast<int> (Ns));
    tPartBuf->SetNumberQ(static_cast<int> (Nas));
    tPartBuf->SetNumberQ(static_cast<int> (Nc));
    tPartBuf->SetNumberQ(static_cast<int> (Nac));
    tPartBuf->SetPDGCode(static_cast<int> (MC));
    aDB->AddParticleType(tPartBuf);
    delete iss;
  }
