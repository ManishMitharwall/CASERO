def Print_input_parameters(Nst,Nel,Sz,max_oc,Processors):
    if type(Sz) == int:
        Sz = Sz % 2
    print('\n*******************INPUT  PARAMETERS*************************')
    print(f'Number of processor used is                     {Processors}')
    print(f'Maximum occupation of each site is              {max_oc}')
    print(f'Number of sites are                             {Nst}')
    print(f'Total number of electrons are                   {Nel}')
    print(f'Sz value is                                     {Sz}')
    print('*************************************************************')

def Print_energies_and_S_value(E,sval):
    print('\n\tEnergy(ev)\t Spin')
    for i in range(len(E)):
        print('\t%.4f \t %.2f ' % (E[i],sval[i]))
    print('Energy gap of ground state and excited state is     %.4f ' %(E[1]-E[0]))

def Print_time_elepsed(st,et):
    hours, rem = divmod(et-st, 3600)
    minutes, seconds = divmod(rem, 60)
    print("\nTime Taken: {:0>2}H:{:0>2}M:{:05.2f}S\n".format(int(hours),int(minutes),seconds))

def Print_Title():
    title = '''
+---------------------------------------------------------------------------------+
|                              D S P M M                                          |
|                           ================                                      |
|                        A      CASCI      CODE                                   |
+---------------------------------------------------------------------------------+'''
    print(title)

