pip install pyiron
create ~/.pyiron


%pip install pyironbase==0.0.3
%Requirement already satisfied: pyironbase==0.0.3 in /Users/glensk/miniconda3/lib/python3.6/site-packages
%Requirement already satisfied: matplotlib in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: tables in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: sqlalchemy in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: h5py in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: pandas in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: six in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: numpy in /Users/glensk/miniconda3/lib/python3.6/site-packages (from pyironbase==0.0.3)
%Requirement already satisfied: python-dateutil>=2.1 in /Users/glensk/miniconda3/lib/python3.6/site-packages (from matplotlib->pyironbase==0.0.3)
%Requirement already satisfied: kiwisolver>=1.0.1 in /Users/glensk/miniconda3/lib/python3.6/site-packages (from matplotlib->pyironbase==0.0.3)
%Requirement already satisfied: cycler>=0.10 in /Users/glensk/miniconda3/lib/python3.6/site-packages (from matplotlib->pyironbase==0.0.3)
%Requirement already satisfied: pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 in /Users/glensk/miniconda3/lib/python3.6/site-packages (from matplotlib->pyironbase==0.0.3)
%Requirement already satisfied: pytz in /Users/glensk/miniconda3/lib/python3.6/site-packages (from matplotlib->pyironbase==0.0.3)
%Requirement already satisfied: numexpr>=2.5.2 in /Users/glensk/miniconda3/lib/python3.6/site-packages (from tables->pyironbase==0.0.3)
%Requirement already satisfied: setuptools in /Users/glensk/miniconda3/lib/python3.6/site-packages (from kiwisolver>=1.0.1->matplotlib->pyironbase==0.0.3)
%
%
DANACH WAR matplotlib kaputt, habe es repairiert mit

conda install matplotlib

