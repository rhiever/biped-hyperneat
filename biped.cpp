/*

Biped HyperNEAT project

Copyright 2013 Randal S. Olson, Joel Lehman, Kenneth Stanley.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

*/

// Test for cylinder vs sphere, by Bram Stolk

#include <ode/odeconfig.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <unistd.h>


#include "rtneat/neat.h"
#include "rtneat/organism.h"
#include "rtneat/noveltyset.h"
#include "rtneat/datarec.h"

#include "simplega/ga.h"
#include "simplega/ConfigFile.h"

static char startgenes_fn[100]="bipedstartgenes";
static char substrate_fn[100] = "substrategenes";
static int novelty_function = 0;

//#define RIGIDITY 

#define NF_COG 14
#define NF_COGSAMP 15
#define NF_LEG 16
#define NF_LEGSAMP 17
#define NF_SS 18

#define NF_FITCU 0
#define NF_FITCUSAMP 1
#define NF_COGCU 2
#define NF_COGCUSAMP 3
#define NF_LEGCU 4
#define NF_LEGCUSAMP 5
#define NF_COGSQ 6
#define NF_COGSAMPSQ 7
#define NF_LEGSQ 8
#define NF_LEGSAMPSQ 9
#define NF_FITSQ 10
#define NF_FITSQSAMP 11
#define NF_FIT 12
#define NF_FITSAMP 13

using namespace std;

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <ode/ode.h>

#ifdef _MSC_VER
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif

#ifdef GRAPHICS
#include <drawstuff/drawstuff.h>
#include "texturepath.h"
#endif

class WalkerDomain;
static WalkerDomain *domain;
static SGA::Population *pop;
static SGA::GA *ga;
static dReal params[]={ 0.2, 3.0, 0.1,     0.15, 0.0, 0.1 ,   0.2,0.3,0.35};


class SineController;
class DummyController;
class Controller;

static Controller* controller;

//NEAT + NS stuff
inline float dist(float x1, float y1, float x2, float y2)
{
	float xd = x1-x2;
	float yd = y1-y2;
	return xd*xd+yd*yd;
}

static	void calculate_delta(dVector3 v1, dVector3 v2, dVector3 o)
	{
		for(int x=0;x<3;x++)
			o[x]=v2[x]-v1[x];
	}

static	void calculate_power(dVector3 v,int pow)
	{
		for(int x=0;x<3;x++)
		{
			float temp=v[x];
			bool sign=false;
			if(temp<0.0)
				sign=true;
			for(int k=1;k<pow;k++)
				v[x]*=temp;
			if(sign)
				v[x]=(-v[x]);
		}
	}
float feet_distance2(vector<float>& t1, vector<float>& x1, vector<float>& y1, vector<float>& t2, vector<float> &x2, vector<float> &y2)
{
	float td=0.0;
	float dummy_distance=0.05;
	int c1=0;
	int c2=0;
	int size1 = x1.size();
	int size2 = x2.size();
	int loopsize=0;
	int bigsize=0;
	float  step_factor=1.0;
	float inc_factor=1.0;

	//cout << x1.size() << " " << y1.size() << " " << x2.size() << " " << y2.size() << endl;
	if(size1<size2)
	{
		loopsize=size1;
		bigsize=size2;	
	}
	else
	{
		loopsize=size2;
		bigsize=size1;
	}
	for(int x=0;x<loopsize;x++)
	{
	//	cout << x << " " <<  y1.size() << " " << y2.size() << endl;
		td+=step_factor*dist(x1[x],y1[x],x2[x],y2[x]);
		//step_factor*=inc_factor;
	}

	for(int x=loopsize;x<bigsize;x++)
	{
	//	cout << x << " " << loopsize << " " << bigsize << endl;
		if(size1<size2)
			td+=step_factor*dist(x2[x],y2[x],0.0,0.0);
		else
			td+=step_factor*dist(x1[x],y1[x],0.0,0.0);
		//step_factor*=inc_factor;
	}

	return td;
}

float feet_distance(vector<float>& t1, vector<float>& x1, vector<float>& y1, vector<float>& t2, vector<float> &x2, vector<float> &y2)
{
	float td=0.0;
	int c1=0;
	int c2=0;
	int size1 = x1.size();
	int size2 = x2.size();

	int time=0;

	float  step_factor=1.0;
    	int step_cnt=0;
    	float dummy_distance=0.05;
    	float inc_factor=1.00;

	while(c1<size1 && c2<size2)
	{
		if(t1[c1+1]<t2[c2+1])
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t1[c1+1]-time);
			time=t1[c1+1];
			c1++;
		}
		else if (t2[c2+1]<t1[c1+1])
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t2[c2+1]-time);
			time=t2[c2+1];
			c2++;
		}
		else	
		{
			td+=step_factor*dist(x1[c1],y1[c1],x2[c2],y2[c2])*(t2[c2+1]-time);
			time=t2[c2+1];
			c1++;
			c2++;
		}
		if(c1 > step_cnt)
		{
			step_cnt=c1;
			step_factor*=inc_factor;
		}
		else if(c2>step_cnt)
		{
			step_cnt=c2;
			step_factor*=inc_factor;
		}
	}

	vector<float>* t;
	vector<float>* x;
	vector<float>* y;
	int c;

	if (t1[t1.size()-1] < t2[t2.size()-1])
	{
		t=&t2;
		x=&x2;
		y=&y2;
		c=c2;
	}
	else
	{
		t=&t1;
		x=&x1;
		y=&y1;
		c=c1;
	}

	step_cnt=c;
	step_factor=pow(inc_factor,c);

	for(int k=0;k<t->size()-1;k++)
	{
		td+=step_factor*dummy_distance*((*t)[k+1]-time);
		time=(*t)[k+1];
		step_factor*=inc_factor;
	}

	if(td<0.0)
	{
		cout << "Wtf" << endl;
		return 0.001;
	}
	return td;	
}

//novelty metric for maze simulation
float walker_novelty_metric(noveltyitem* x,noveltyitem* y)
{
	float dist=0.0;
	if(novelty_function!=NF_SS)
	{

	int size = x->data[0].size();
	
	for(int k=0;k<size;k++)
	{
		float delta = x->data[0][k]-y->data[0][k];
		dist+=delta*delta;
	}

	return dist;	

	}
	else
	{
	float left_feet;
	float right_feet;
	left_feet = feet_distance2(x->data[0],x->data[1],x->data[2],y->data[0],y->data[1],y->data[2]);
	right_feet = feet_distance2(x->data[3],x->data[4],x->data[5],y->data[3],y->data[4],y->data[5]);
	
	return left_feet+right_feet;
	}
}

void biped_epoch(NEAT::Population *pop,bool novelty=false);
void biped_realtime_loop(NEAT::Population *pop,bool novelty=false);

NEAT::Population *biped_realtime(bool novelty=false);
noveltyitem* biped_evaluate(NEAT::Organism* org,data_record* data=NULL);

//globals
static NEAT::Population *neatpop;
static noveltyarchive archive(1.0,*walker_novelty_metric);
static data_rec Record; //stores run information
static int indiv_counter=0;
static bool Novelty=false;
static char outdir[50]="";

class Controller
{
public:
	vector<dReal> outs;
	int size;
	bool scale;
	bool debug;
	ofstream* dbgfile;
	Controller() { } 
	Controller(int s,bool deb=false)
	{
		debug=deb;
		scale=true;
		size=s;
		for(int x=0;x<size;x++)
			outs.push_back(0.0);
		if (debug)
			dbgfile=new ofstream("controller.log");
	}
	virtual void update(double time,vector<dReal> sensors) {
		if(debug)
                {
		for(int x=0;x<size;x++)
			{
			 *dbgfile << outs[x] << " ";
			}
			*dbgfile << "\n";
		}
	}
	virtual vector<dReal>* get_outputs()
	{
		return &outs;
	}
	virtual ~Controller()
	{
		if(debug)
			delete dbgfile;
	}
};

class DummyController: public Controller
{
public:
	DummyController():Controller(6,false)
	{
		scale=false;
		for(int x=0;x<size;x++)
			outs[x]=0.0;
	}
};

class CTRNNController: public Controller
{
public:
	NEAT::Network* net;
	CTRNNController(NEAT::Network* n,bool dbg=false):Controller(6,dbg)
	{
		scale=true;
		net=n;
		n->init_ctrnn();
		
		for(int x=0;x<50;x++)
			net->activate_ctrnn(0.02);
		for(int x=0;x<size;x++)
			outs[x]=net->outputs[x]->output;
	}
	virtual void update(double time,vector<dReal> sensors) 
	{
		double sens[20];
		for(int x=0;x<sensors.size();x++)
			sens[x]=sensors[x];
		net->load_sensors(sens);
		net->activate_ctrnn(time);
		for(int x=0;x<6;x++)
		{
			
			outs[x]=net->outputs[x]->output;
		}
		Controller::update(time,sensors);
	}
};

class SineController: public Controller
{
public:
	double cum_time;
	vector<dReal> params;
	SineController(int size):Controller(size)
	{
		scale=false;
		cum_time=0.0;
	}
	void load_params(dReal *p)
	{
		for(int x=0;x<(size/2)*3;x++)
		{
			params.push_back(p[x]);
		}
	}
	virtual void update(double time,vector<dReal> sensors) 
	{
		for(int x=0;x<size;x+=2)
		{
			double constant=params[3*(x/2)];
			double phase=params[3*(x/2)+1];
			double amplitude=params[3*(x/2)+2];
			outs[x]=constant+amplitude*sin(3.14*2.0*cum_time+phase);
			outs[x+1]=constant+amplitude*sin(3.14*2.0*cum_time+phase+3.14);
		}
		cum_time+=time;
	}
};

class Creature
{
public:
	bool movie_rec;
	bool movie_play;

	ofstream* movie;
	ofstream* movie_rot;

	ifstream* movie_in;
	ifstream* movie_rot_in;

	Controller* controller;
	vector<dGeomID> geoms;
	vector<dBodyID> bodies;
	vector<bool> onground;
	vector<dJointID> joints;
	
	vector<dReal> current_angles;
	vector<dReal> desired_angles;
	
	vector<dReal> lo_limit;
	vector<dReal> hi_limit;

	vector<dReal> sensors;
	vector<dReal> desired_angvel;
	vector<dReal> delta_angles;

	vector<dReal> p_terms;
	vector<dReal> d_terms;

	dWorldID world;
	dSpaceID space;
	dVector3 pos;
	
	Creature(bool mov=false,bool play=false) {
	if(play)
		mov=false;
	movie_play=play;
	movie_rec=mov;
	if(play)
	{
		movie_in=new ifstream("movie_pos.dat");
		movie_rot_in=new ifstream("movie_rot.dat");
	}
	if(mov)
	{
		movie=new ofstream("movie_pos.dat");
		movie_rot=new ofstream("movie_rot.dat");
	}

	}
	
	virtual void Update(double timestep)
	{
		if(movie_rec)
		{

		for(int x=0;x<geoms.size();x++)
		{
			dQuaternion rot;
			const dReal* pos=dGeomGetPosition(geoms[x]);
			dGeomGetQuaternion(geoms[x],rot);
			*movie << pos[0] << " " << pos[1] << " " << pos[2] << " ";
			*movie_rot << rot[0] << " " << rot[1] << " " << rot[2] << " " << rot[3] << " ";
		}
		*movie << endl;
		*movie_rot << endl;
		}
		if(movie_play)
		{
		for(int k=0;k<geoms.size();k++)
		{
		 dReal x,y,z;
		 dQuaternion q;
		 *movie_in >> x >> y >> z;
		 *movie_rot_in >> q[0] >> q[1] >> q[2] >> q[3];
		 dGeomSetPosition(geoms[k],x,y,z);
		 dGeomSetQuaternion(geoms[k],q);
		}

		}
	}

	dReal TotalMass()
	{
		dReal total_mass=0.0;
		for(int x=0;x<bodies.size();x++)
		{
			dMass m;
			dBodyGetMass(bodies[x],&m);
			total_mass+=m.mass;
		}
		return total_mass;
	}

	void CenterOfMass(dVector3 center)
	{
		dReal total_mass=0.0;
		dVector3 accum={0.0,0.0,0.0};
		const dReal* bpos;
		for(int x=0;x<bodies.size();x++)
		{
			dMass m;
			dBodyGetMass(bodies[x],&m);
			bpos=dBodyGetPosition(bodies[x]);
			total_mass+=m.mass;
			for(int y=0;y<3;y++)
				accum[y]+=m.mass*bpos[y];
		}

		for(int x=0;x<3;x++)
			center[x]=accum[x]/total_mass;
	}

	virtual void Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont)
	{
		controller=cont;
		world=worldi;
		space=spacei;
		pos[0]=posi[0];
		pos[1]=posi[1];
		pos[2]=posi[2];
	}

	virtual void Destroy()
	{
		for(int x=0;x<geoms.size();x++)
		  dGeomDestroy(geoms[x]);
	}

	virtual bool abort()=0;
	virtual dReal fitness()=0;
	int add_fixed(int b1,int b2)
	{
		dBodyID bd1=bodies[b1];
		dBodyID bd2;
		if (b2!=(-1))
		{
			bd2=bodies[b2];
		}
		else
		{
			bd2=0;
		}
		dJointID tempjoint = dJointCreateFixed(world,0);
		dJointAttach(tempjoint,bd1,bd2);
		dJointSetFixed(tempjoint);
		joints.push_back(tempjoint);
		return joints.size()-1;
	}
	int add_hinge(int b1,int b2,dVector3 anchor,dVector3 axis,dReal lostop, dReal histop, dReal fmax)
	{
		dJointID tempjoint = dJointCreateHinge(world,0);
		dJointAttach(tempjoint,bodies[b1],bodies[b2]);
		dJointSetHingeAnchor(tempjoint,pos[0]+anchor[0],pos[1]+anchor[1],pos[2]+anchor[2]);
		dJointSetHingeAxis(tempjoint,axis[0],axis[1],axis[2]);
		dJointSetHingeParam(tempjoint,dParamLoStop,lostop);
		dJointSetHingeParam(tempjoint,dParamHiStop,histop);
		dJointSetHingeParam(tempjoint,dParamFMax,fmax);
		joints.push_back(tempjoint);
		
		return joints.size()-1;
	}

	int add_universal(int b1,int b2,dVector3 anchor, dVector3 axis1, dVector3 axis2, dReal lostop1, dReal histop1, dReal lostop2, dReal histop2, dReal fmax1, dReal fmax2)
	{
		dJointID tempjoint = dJointCreateUniversal(world,0);
		dJointAttach(tempjoint,bodies[b1],bodies[b2]);
		dJointSetUniversalAnchor(tempjoint,pos[0]+anchor[0],pos[1]+anchor[1],pos[2]+anchor[2]);
		dJointSetUniversalAxis1(tempjoint,axis1[0],axis1[1],axis1[2]);
    	dJointSetUniversalAxis2(tempjoint,axis2[0],axis2[1],axis2[2]);
		dJointSetUniversalParam(tempjoint,dParamLoStop,lostop1);
		dJointSetUniversalParam(tempjoint,dParamHiStop,histop1);
		dJointSetUniversalParam(tempjoint,dParamLoStop2,lostop2);
		dJointSetUniversalParam(tempjoint,dParamHiStop2,histop2);
		dJointSetUniversalParam(tempjoint,dParamFMax,fmax1);
		
		
  	    dJointSetUniversalParam(tempjoint,dParamFMax2,fmax2);

		joints.push_back(tempjoint);
		return joints.size()-1;
	}
	int add_box(dReal density, dReal lx, dReal ly, dReal lz, const dVector3 p)
	{
		dBodyID tempbody;
		dGeomID tempgeom;
		dMass m;
   	    tempbody = dBodyCreate (world);
		dMassSetBoxTotal(&m,density,lx,ly,lz);
		tempgeom = dCreateBox(0,lx,ly,lz);
		dGeomSetBody(tempgeom,tempbody);
		dBodySetPosition(tempbody,pos[0]+p[0],pos[1]+p[1],pos[2]+p[2]);
		dSpaceAdd(space,tempgeom);

		 bodies.push_back(tempbody);
		 geoms.push_back(tempgeom);
		 onground.push_back(false);
		 return bodies.size()-1;
	}

	int add_sphere(dReal density, dReal radius, const dVector3 p)
	{
		
  		 dBodyID tempbody;
		 dGeomID tempgeom;
		 dMass m;
		 tempbody = dBodyCreate (world);
		 dMassSetSphereTotal(&m,density,radius);
		 tempgeom = dCreateSphere(0,radius);
		 dGeomSetBody(tempgeom,tempbody);
		 dBodySetMass(tempbody,&m);
		 dBodySetPosition(tempbody,pos[0]+p[0],pos[1]+p[1],pos[2]+p[2]);
		 dSpaceAdd(space,tempgeom);

		 bodies.push_back(tempbody);
		 geoms.push_back(tempgeom);
		 onground.push_back(false);
		 return bodies.size()-1;
		
	}

	int add_cylinder(int axis, dReal density,dReal length, dReal radius, const dVector3 p,dBodyID* k=NULL)
	{
		dReal a[]={0.0,0.0,0.0};

		 dBodyID tempbody;
		 dGeomID tempgeom;
		 dMass m;

		tempbody = dBodyCreate (world);
		if(k!=NULL)
			(*k)=tempbody;
		dQuaternion q;
		if (axis==1)
		{
			a[1]=1.0;
		}
		else if (axis==2)
		{
			a[0]=1.0;
		}
		else
		{
			a[2]=1.0;
		}
	    dQFromAxisAndAngle (q,a[0],a[1],a[2], M_PI * 0.5);
	    dBodySetQuaternion (tempbody,q);
		dMassSetCylinderTotal (&m,density,axis,radius,length);
		dBodySetMass (tempbody,&m);
		tempgeom = dCreateCylinder(0, radius, length);
		dGeomSetBody (tempgeom,tempbody);
		dBodySetPosition (tempbody, pos[0]+p[0],pos[1]+p[1], pos[2]+p[2]);		
		dSpaceAdd (space, tempgeom);

		geoms.push_back(tempgeom);
		bodies.push_back(tempbody);
		onground.push_back(false);
		return bodies.size()-1;
	}
	virtual ~Creature()
	{

	if(movie_rec)
	{
		delete movie;
		delete movie_rot;
		cout << "terminating.." << endl;
	}
	if(movie_play)
	{
		delete movie_in;
		delete movie_rot_in;
	}

	}
};

//BIPED PARAMETERS
static dReal SCALE_FACTOR = 3;
static dReal FOOTX_SZ =1.0/SCALE_FACTOR;
static dReal FOOTY_SZ =0.5/SCALE_FACTOR;
static dReal FOOTZ_SZ =1.0/SCALE_FACTOR;
static dReal LLEG_LEN =1.0/SCALE_FACTOR;
static dReal LLEG_RAD =0.2/SCALE_FACTOR;
static dReal ULEG_LEN =1.0/SCALE_FACTOR;
static dReal ULEG_RAD =0.2/SCALE_FACTOR;
static dReal TORSO_LEN =1.0/SCALE_FACTOR;
static dReal TORSO_RAD =0.3/SCALE_FACTOR;
static dReal ORIG_HEIGHT= (TORSO_RAD/2.0+ULEG_LEN+LLEG_LEN+FOOTZ_SZ);
static dReal DENSITY=0.5;
static dReal TORSO_DENSITY=1.0;
static dReal FOOT_DENSITY=0.1;
static dReal MAXTORQUE_FOOT= 10.0;
static dReal MAXTORQUE_KNEE= 5.0;
static dReal MAXTORQUE_HIPMINOR= 5.0;
static dReal MAXTORQUE_HIPMAJOR= 5.0;
static dReal P_CONSTANT= 9.0;
static dReal D_CONSTANT= 0.0;
static dReal FOOTFACTOR= 5.0;

static bool bDisplay = false;
static bool bDoEvolution = true;
static bool bSlowDown = false;
#ifdef GRAPHICS
static bool bNoVis=true;//false;
#else
static bool bNoVis=true;
#endif
static bool bMoviePlay=false;

void load_parameters()
{
	ConfigFile config( "config.txt" );
	DENSITY=config.read<dReal>("DENSITY");
	P_CONSTANT=config.read<dReal>("P_CONSTANT");
	D_CONSTANT=config.read<dReal>("D_CONSTANT");
	MAXTORQUE_KNEE=config.read<dReal>("MAXTORQUE_KNEE");
	MAXTORQUE_HIPMINOR=config.read<dReal>("MAXTORQUE_HIPMINOR");
	MAXTORQUE_HIPMAJOR=config.read<dReal>("MAXTORQUE_HIPMAJOR");
	MAXTORQUE_FOOT=config.read<dReal>("MAXTORQUE_FOOT");
	//FOOT_P_CONSTANT=config.read<dReal>("FOOT_P_CONSTANT");
	//FOOT_D_CONSTANT=config.read<dReal>("FOOT_D_CONSTANT");
	int novis=config.read<int>("NOVIS");
	if(novis==0)
		bNoVis=false;
	else
		bNoVis=true;
	bMoviePlay=(bool)config.read<int>("MOVIE");
}

dReal scaleval(dReal v, dReal min, dReal max)
{
if(v<min)
	return 0.0;
if(v>max)
	return 1.0;
return (v-min)/(max-min);
}

class Biped: public Creature
{
public:
	int step;
	dVector3 orig_com;
	dVector3 orig_left;
	dVector3 orig_right;
	dVector3 curr_com;
	bool log;
	ofstream* logfile;
	dJointFeedback feedback[6];

	//keeping track of foot positionxorz
	bool leftdown;
	bool rightdown;
	
	bool leftrigid;
	bool rightrigid;

	int lastdown; //which foot was last down

	vector<float> lft; //left foot time
	vector<float> lfx; //x
	vector<float> lfy; //y
	vector<float> rft; //right foot time
	vector<float> rfx; //x
	vector<float> rfy; //y

 Biped(bool logging=false,bool movie=false):Creature(logging,movie) {
	step=0;
	
	leftdown=false;
	rightdown=false;
	
	leftrigid=false;
	rightrigid=false;

	lastdown=0;

	log=logging;

	 if(log)
	 {
		 logfile=new ofstream("log.dat");
	 }
	 
	for(int x=0;x<6;x++)
	{
		p_terms.push_back(P_CONSTANT);
		d_terms.push_back(D_CONSTANT);
		desired_angles.push_back(0.0);
		current_angles.push_back(0.0);
		delta_angles.push_back(0.0);
		desired_angvel.push_back(0.0);
		lo_limit.push_back(0.0);
		hi_limit.push_back(0.0);
	}

	for(int x=0;x<8;x++)
		sensors.push_back(0.0);
 } 
	int add_foot(dReal density, dReal radius, const dVector3 p)
	{
		
  		 dBodyID tempbody;
		 dGeomID maingeom[3];
		 dGeomID tempgeom[3];

		 dMass m;
		 tempbody = dBodyCreate (world);
		 dMassSetSphereTotal(&m,density,radius);
		 dBodySetPosition(tempbody,pos[0]+p[0],pos[1]+p[1],pos[2]+p[2]);
		 bodies.push_back(tempbody);
		 onground.push_back(false);
		 
		 for(int x=0;x<3;x++)
		 {
                 maingeom[x] = dCreateGeomTransform(space);
		 tempgeom[x] = dCreateSphere(0,radius);
		 
		 dGeomSetBody(maingeom[x],tempbody);
		 
		 dGeomTransformSetGeom(maingeom[x],tempgeom[x]);
		 dGeomTransformSetCleanup(maingeom[x],1);
		 dGeomTransformSetInfo(maingeom[x],1);
		 
		 geoms.push_back(maingeom[x]);
		 
		 if(x==0)
		  dGeomSetPosition(tempgeom[x],0.07,-0.08,0.0);
		 if(x==1)
                  dGeomSetPosition(tempgeom[x],0.07,0.08,0.0);
		 if(x==2)
                  dGeomSetPosition(tempgeom[x],-0.20,0.0,0.0);

		 }
		 return bodies.size()-1;
	}	
		
		
 virtual dReal fitness()
 {
	 dVector3 new_com;
	 CenterOfMass(new_com);
	 double fitness=0.0;
	 for(int x=0;x<2;x++)
	 {
		 double delta=new_com[x]-orig_com[x];
		 delta*=delta;
		 fitness+=delta;
	 }
	 return sqrt(fitness);
 }
	
	~Biped()
	{
		if(log)
			delete logfile;
	}

 virtual bool abort()
 {
	 const dReal* torsoPos;
	 torsoPos=dBodyGetPosition(bodies[6]);
	 if (torsoPos[2]< 0.5*(ORIG_HEIGHT))
		 return true;
	 return false;
 }

 void create_leg(dVector3 offset)
 {
	 dVector3 xAxis={1.0,0.0,0.0};
	 dVector3 yAxis={0.0,-1.0,0.0};
	 dVector3 zAxis={0.0,0.0,1.0};

	 dVector3 p={offset[0],offset[1],offset[2]};

	 //dVector3 foot_pos = {p[0]-0.3*FOOTX_SZ,p[1],p[2]+(FOOTZ_SZ/2.0)};
	 //int foot=add_box(DENSITY,FOOTX_SZ,FOOTY_SZ,FOOTZ_SZ,foot_pos);

	 dVector3 foot_pos = {p[0],p[1],p[2]+(FOOTZ_SZ/2.0)};
	 
	 int foot=add_sphere(FOOT_DENSITY,FOOTZ_SZ/2.0,foot_pos);
	
	 
	 dVector3 lower_pos = {p[0],p[1],p[2]+FOOTZ_SZ+LLEG_LEN/2.0};
	 int lowerleg = add_cylinder(3,DENSITY,LLEG_LEN,LLEG_RAD,lower_pos);
	 dVector3 upper_pos = {p[0],p[1],p[2]+FOOTZ_SZ+LLEG_LEN+ULEG_LEN/2.0};
	 int upperleg = add_cylinder(3,DENSITY,ULEG_LEN,ULEG_RAD,upper_pos);

	 dVector3 foot_joint_a = {p[0],p[1],p[2]+(FOOTZ_SZ)};
	 dVector3 knee_joint_a = {p[0],p[1],p[2]+FOOTZ_SZ+LLEG_LEN};
	 
	 add_fixed(foot,lowerleg);
	 //add_universal(foot,lowerleg,foot_joint_a,xAxis,yAxis,-0.1,0.1,-0.1,0.1,MAXTORQUE_FOOT,MAXTORQUE_FOOT);
	 add_hinge(lowerleg,upperleg,knee_joint_a,yAxis,-1.4,0.0,MAXTORQUE_KNEE);
 }

 virtual void Create(dWorldID worldi, dSpaceID spacei, dVector3 posi,Controller* cont)
 {
	 dVector3 xAxis={1.0,0.0,0.0};
	 dVector3 nxAxis={-1.0,0.0,0.0};
	 dVector3 yAxis={0.0,1.0,0.0};
	 dVector3 zAxis={0.0,0.0,1.0};
	
	 Creature::Create(worldi,spacei,posi,cont);
	 dVector3 leftLegPos={0.0,0.0,0.0};
	 dVector3 rightLegPos={0.0,TORSO_LEN+ULEG_RAD,0.0};
	 //dVector3 torsoPos={0.0,(TORSO_LEN+ULEG_RAD)/2.0,ULEG_LEN+LLEG_LEN+FOOTZ_SZ};
	 dVector3 torsoPos={0.0,(TORSO_LEN+ULEG_RAD)/2.0,ULEG_LEN+LLEG_LEN+FOOTZ_SZ};
	

	 dVector3 leftHip={leftLegPos[0],leftLegPos[1]+ULEG_RAD,torsoPos[2]};
	 dVector3 rightHip={rightLegPos[0],rightLegPos[1]-ULEG_RAD,torsoPos[2]};

	 create_leg(leftLegPos);
	 create_leg(rightLegPos);

	 int torso=add_cylinder(2,TORSO_DENSITY,TORSO_LEN,TORSO_RAD,torsoPos);
							//-1.4,1.4
	 add_universal(torso,2,leftHip,xAxis,yAxis,-0.8,0.8,-1.3,1.6,MAXTORQUE_HIPMINOR,MAXTORQUE_HIPMAJOR);
	 add_universal(torso,5,rightHip,nxAxis,yAxis,-0.8,0.8,-1.3,1.6,MAXTORQUE_HIPMINOR,MAXTORQUE_HIPMAJOR);
	 
	 lo_limit[0]=dJointGetHingeParam(joints[1],dParamLoStop);
	 lo_limit[1]=dJointGetHingeParam(joints[3],dParamLoStop);	 
	 hi_limit[0]=dJointGetHingeParam(joints[1],dParamHiStop);
	 hi_limit[1]=dJointGetHingeParam(joints[3],dParamHiStop);

	 lo_limit[2]=dJointGetUniversalParam(joints[4],dParamLoStop);
	 lo_limit[3]=dJointGetUniversalParam(joints[5],dParamLoStop);
	 hi_limit[2]=dJointGetUniversalParam(joints[4],dParamHiStop);
	 hi_limit[3]=dJointGetUniversalParam(joints[5],dParamHiStop);

	 lo_limit[4]=dJointGetUniversalParam(joints[4],dParamLoStop2);
	 lo_limit[5]=dJointGetUniversalParam(joints[5],dParamLoStop2);
	 hi_limit[4]=dJointGetUniversalParam(joints[4],dParamHiStop2);
	 hi_limit[5]=dJointGetUniversalParam(joints[5],dParamHiStop2);

/*
	 lo_limit[6]=dJointGetUniversalParam(joints[0],dParamLoStop);
	 lo_limit[7]=dJointGetUniversalParam(joints[2],dParamLoStop);
	 hi_limit[6]=dJointGetUniversalParam(joints[0],dParamHiStop);
	 hi_limit[7]=dJointGetUniversalParam(joints[2],dParamHiStop);

	 lo_limit[8]=dJointGetUniversalParam(joints[0],dParamLoStop2);
	 lo_limit[9]=dJointGetUniversalParam(joints[2],dParamLoStop2);
	 hi_limit[8]=dJointGetUniversalParam(joints[0],dParamHiStop2);
	 hi_limit[9]=dJointGetUniversalParam(joints[2],dParamHiStop2);
*/
	 dJointSetFeedback(joints[1],&feedback[0]);
	 dJointSetFeedback(joints[3],&feedback[1]);
	 dJointSetFeedback(joints[4],&feedback[2]);
	 dJointSetFeedback(joints[5],&feedback[3]);
	 dJointSetFeedback(joints[0],&feedback[4]);
	 dJointSetFeedback(joints[2],&feedback[5]);

	 //add_fixed(torso,-1);
		 
	CenterOfMass(orig_com);
	CenterOfMass(curr_com);
	orig_left[0]=dBodyGetPosition(bodies[0])[0];
	orig_right[0]=dBodyGetPosition(bodies[0])[1];	
	orig_left[1]=dBodyGetPosition(bodies[3])[0];
	orig_right[1]=dBodyGetPosition(bodies[3])[1];
	orig_left[2]=0.0;
	orig_right[2]=0.0;
 }

 void print_behavior()
 {

	cout << "LEFTFOOTSTEPS: " << lft.size() << endl;
	cout << "RIGHTFOOTSTEPS: " << rft.size() << endl;

	for(int x=0;x<lft.size();x++)
	{
		cout << "LFT " << x << " time: " << lft[x] << " x: " << lfx[x] << " y: " << lfy[x] << endl;
	}
	for(int x=0;x<rft.size();x++)
	{
		cout << "RFT " << x << " time: " << rft[x] << " x: " << rfx[x] << " y: " << rfy[x] << endl;
	}

 }

 // 2 input
 virtual void Update(double timestep)
 {
	Creature::Update(timestep);
	if (movie_play)
		return;
	dReal old_angles[10];

	step++;

	for(int x=0;x<6;x++)
		 old_angles[x]=current_angles[x];
	 
	 //read current angles
	current_angles[0]=dJointGetHingeAngle(joints[1]); //left knee
	current_angles[1]=dJointGetHingeAngle(joints[3]); //right knee
	
	current_angles[2]=dJointGetUniversalAngle1(joints[4]); //left outhip
	current_angles[3]=dJointGetUniversalAngle1(joints[5]); //right outhip

	current_angles[4]=dJointGetUniversalAngle2(joints[4]); //left mainhip
	current_angles[5]=dJointGetUniversalAngle2(joints[5]); //right mainhip


	for(int x=0;x<6;x++)
		 delta_angles[x]=(current_angles[x]-old_angles[x])/timestep;

	//record behavior
	bool newleftdown=onground[0];
	bool newrightdown=onground[3];	 

	const dReal* q = dBodyGetQuaternion(bodies[6]);
	dReal tanyaw = 2.0*(q[0]*q[1]+q[3]*q[2])/
			(q[3]*q[3]+q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
	dReal sinpitch = -2.0*(q[0]*q[2]-q[3]*q[1]);

	dReal tanroll = 2.0*(q[3]*q[0]+q[1]*q[2])/
			(q[3]*q[3]-q[0]*q[0]-q[1]*q[1]+q[2]*q[2]);

	dReal yaw=atan(tanyaw);
	dReal pitch=asin(sinpitch);
	dReal roll=atan(tanroll);

	//angle sensors
	sensors[2]=scaleval(yaw,-3.14/2,3.14/2);
	sensors[3]=scaleval(pitch,-3.14/2,3.14/2);
	sensors[4]=scaleval(roll,-3.14/2,3.14/2);
	
	//accelerometers
	/*
	sensors[2]=dBodyGetAngularVel(bodies[6])[0];
	sensors[3]=dBodyGetAngularVel(bodies[6])[1];
	sensors[4]=dBodyGetAngularVel(bodies[6])[2];
 	*/

	

	for(int x=0;x<6;x++)
		sensors[2+x]=current_angles[x];
	

	if(newleftdown)
		sensors[0]=1.0;
	else
		sensors[0]=0.0;

	if(newrightdown)
		sensors[1]=1.0;
	else
		sensors[1]=0.0;

	//update controller
	controller->update(timestep,sensors);
	vector<dReal>* outs=controller->get_outputs();

	//calculate values to pass to joint motors
	if(log)
	 (*logfile) << "-PID" <<endl;
	


	for(int x=0;x<6;x++)
	{
	
 		desired_angles[x]=(*outs)[x];
	
		
		if(controller->scale)
		{
			if(desired_angles[x]>1.0) desired_angles[x]=1.0;
			if(desired_angles[x]<0.0) desired_angles[x]=0.0;
			desired_angles[x]= lo_limit[x]+(hi_limit[x]-lo_limit[x])*desired_angles[x];
		}
	}

	for(int x=0;x<6;x++)
	{
		dReal delta=desired_angles[x]-current_angles[x];


		double p_term = p_terms[x]* delta;
		double d_term = (-d_terms[x]*delta_angles[x]);
		desired_angvel[x]=p_term+d_term;
		if(log)
			(*logfile) << p_term << " " << d_term << " " << desired_angvel[x] 
				<< " " << delta_angles[x] << " " << desired_angles[x] << " " << current_angles[x] << endl;
	}
	
	 if(log)
		(*logfile) << "-FEED" <<endl;
	
	 if(log)
	 for(int x=0;x<4;x++)
	 {
		 for(int k=0;k<3;k++)
			 (*logfile) << feedback[x].f1[k] << " ";
		 for (int k=0;k<3;k++)
			 (*logfile) << feedback[x].f2[k] << " ";
		 (*logfile) << endl;
	 }
	
	//update joint motors
	dJointSetHingeParam(joints[1],dParamVel,desired_angvel[0]); //left knee
	dJointSetHingeParam(joints[3],dParamVel,desired_angvel[1]); //right knee

	dJointSetUniversalParam(joints[4],dParamVel,desired_angvel[2]); //left hipout
	dJointSetUniversalParam(joints[5],dParamVel,desired_angvel[3]); //right hipout

	dJointSetUniversalParam(joints[4],dParamVel2,desired_angvel[4]); //left hipmain
	dJointSetUniversalParam(joints[5],dParamVel2,desired_angvel[5]); //right hipmain


	if(!leftdown && newleftdown)
	{
		if(lft.size()==0 || (step-lft[lft.size()-1] > 100 && lastdown!=1))
		{
                       
			CenterOfMass(curr_com);
			lft.push_back(step);
			lfx.push_back(curr_com[0]);
			lfy.push_back(curr_com[1]);
			//lfx.push_back(dBodyGetPosition(bodies[0])[0]);
			//lfy.push_back(dBodyGetPosition(bodies[0])[1]);

			if (novelty_function % 2 == 1)
				lastdown=(1);
			else
				lastdown=0;  //if this is set to 0, we don't care if feet sequence alternates
		}
	}
	if(!rightdown && newrightdown)
	{
		//enable for separate foot tracking
		
		if(rft.size()==0 || (step-rft[rft.size()-1] > 100 && lastdown!=-1))
		{
		
			CenterOfMass(curr_com);
			rft.push_back(step);
			rfx.push_back(curr_com[0]);
			rfy.push_back(curr_com[1]);
			//rfx.push_back(dBodyGetPosition(bodies[3])[0]);
			//rfy.push_back(dBodyGetPosition(bodies[3])[1]);
			if (novelty_function % 2 ==1)
				lastdown=(-1);
			else		
				lastdown=0;  //if this is set to 0, we don't care if feet sequence alternates
		}
		/*
		if(step==1 || (step-lft[lft.size()-1] > 200 && lastdown!=(-1)))
		{
			lft.push_back(step);
			lfx.push_back(dBodyGetPosition(bodies[3])[0]);
			lfy.push_back(dBodyGetPosition(bodies[3])[1]);
			lastdown=(-1);
		}
		*/
	}
	//don't let first recorded instance of both feet down set the lastdown criteria
	if(step==1)
		lastdown=0;

	leftdown=newleftdown;
	rightdown=newrightdown;
	//reset ground sensors for feetz
	onground[0]=false;
	onground[3]=false;

 }

 // 6 input
 virtual void Update_dnu(double timestep)
 {
	Creature::Update(timestep);
	if (movie_play)
		return;
	dReal old_angles[10];

	step++;

	for(int x=0;x<6;x++)
		 old_angles[x]=current_angles[x];
	 
	//read current angles
	/*current_angles[0]=dJointGetHingeAngle(joints[1]); //left knee
	current_angles[1]=dJointGetHingeAngle(joints[3]); //right knee
	
	current_angles[2]=dJointGetUniversalAngle1(joints[4]); //left outhip
	current_angles[3]=dJointGetUniversalAngle1(joints[5]); //right outhip

	current_angles[4]=dJointGetUniversalAngle2(joints[4]); //left mainhip
	current_angles[5]=dJointGetUniversalAngle2(joints[5]);*/ //right mainhip

        current_angles[0]=dJointGetHingeAngle(joints[1]); //left knee
        current_angles[5]=dJointGetHingeAngle(joints[3]); //right knee

        current_angles[1]=dJointGetUniversalAngle1(joints[4]); //left outhip
        current_angles[4]=dJointGetUniversalAngle1(joints[5]); //right outhip

        current_angles[2]=dJointGetUniversalAngle2(joints[4]); //left mainhip
        current_angles[3]=dJointGetUniversalAngle2(joints[5]); //right mainhip



	for(int x=0;x<6;x++)
		 delta_angles[x]=(current_angles[x]-old_angles[x])/timestep;

	//record behavior
	bool newleftdown=onground[0];
	bool newrightdown=onground[3];	 

	const dReal* q = dBodyGetQuaternion(bodies[6]);
	dReal tanyaw = 2.0*(q[0]*q[1]+q[3]*q[2])/
			(q[3]*q[3]+q[0]*q[0]-q[1]*q[1]-q[2]*q[2]);
	dReal sinpitch = -2.0*(q[0]*q[2]-q[3]*q[1]);

	dReal tanroll = 2.0*(q[3]*q[0]+q[1]*q[2])/
			(q[3]*q[3]-q[0]*q[0]-q[1]*q[1]+q[2]*q[2]);

	dReal yaw=atan(tanyaw);
	dReal pitch=asin(sinpitch);
	dReal roll=atan(tanroll);

	//angle sensors
	sensors[2]=scaleval(yaw,-3.14/2,3.14/2);
	sensors[3]=scaleval(pitch,-3.14/2,3.14/2);
	sensors[4]=scaleval(roll,-3.14/2,3.14/2);
	
	//accelerometers
	/*
	sensors[2]=dBodyGetAngularVel(bodies[6])[0];
	sensors[3]=dBodyGetAngularVel(bodies[6])[1];
	sensors[4]=dBodyGetAngularVel(bodies[6])[2];
 	*/

	// Revised for HyperNEAT
	for(int x=0;x<6;x++)
		sensors[x]=current_angles[x];

	/*for(int x=0;x<6;x++)
		sensors[2+x]=current_angles[x];
	

	if(newleftdown)
		sensors[0]=1.0;
	else
		sensors[0]=0.0;

	if(newrightdown)
		sensors[1]=1.0;
	else
		sensors[1]=0.0;*/

	//update controller
	controller->update(timestep,sensors);
	vector<dReal>* outs=controller->get_outputs();

	//calculate values to pass to joint motors
	if(log)
	 (*logfile) << "-PID" <<endl;
	
	for(int x=0;x<6;x++)
	{
	
 		desired_angles[x]=(*outs)[x];
	
		
		if(controller->scale)
		{
			if(desired_angles[x]>1.0) desired_angles[x]=1.0;
			if(desired_angles[x]<0.0) desired_angles[x]=0.0;
			desired_angles[x]= lo_limit[x]+(hi_limit[x]-lo_limit[x])*desired_angles[x];
		}
	}

	for(int x=0;x<6;x++)
	{
		dReal delta=desired_angles[x]-current_angles[x];


		double p_term = p_terms[x]* delta;
		double d_term = (-d_terms[x]*delta_angles[x]);
		desired_angvel[x]=p_term+d_term;
		if(log)
			(*logfile) << p_term << " " << d_term << " " << desired_angvel[x] 
				<< " " << delta_angles[x] << " " << desired_angles[x] << " " << current_angles[x] << endl;
	}
	
	 if(log)
		(*logfile) << "-FEED" <<endl;
	
	 if(log)
	 for(int x=0;x<4;x++)
	 {
		 for(int k=0;k<3;k++)
			 (*logfile) << feedback[x].f1[k] << " ";
		 for (int k=0;k<3;k++)
			 (*logfile) << feedback[x].f2[k] << " ";
		 (*logfile) << endl;
	 }
	
	//update joint motors
	dJointSetHingeParam(joints[1],dParamVel,desired_angvel[0]); //left knee
	dJointSetHingeParam(joints[3],dParamVel,desired_angvel[1]); //right knee

	dJointSetUniversalParam(joints[4],dParamVel,desired_angvel[2]); //left hipout
	dJointSetUniversalParam(joints[5],dParamVel,desired_angvel[3]); //right hipout

	dJointSetUniversalParam(joints[4],dParamVel2,desired_angvel[4]); //left hipmain
	dJointSetUniversalParam(joints[5],dParamVel2,desired_angvel[5]); //right hipmain


	if(!leftdown && newleftdown)
	{
		if(lft.size()==0 || (step-lft[lft.size()-1] > 100 && lastdown!=1))
		{
                       
			CenterOfMass(curr_com);
			lft.push_back(step);
			lfx.push_back(curr_com[0]);
			lfy.push_back(curr_com[1]);
			//lfx.push_back(dBodyGetPosition(bodies[0])[0]);
			//lfy.push_back(dBodyGetPosition(bodies[0])[1]);

			if (novelty_function % 2 == 1)
				lastdown=(1);
			else
				lastdown=0;  //if this is set to 0, we don't care if feet sequence alternates
		}
	}
	if(!rightdown && newrightdown)
	{
		//enable for separate foot tracking
		
		if(rft.size()==0 || (step-rft[rft.size()-1] > 100 && lastdown!=-1))
		{
		
			CenterOfMass(curr_com);
			rft.push_back(step);
			rfx.push_back(curr_com[0]);
			rfy.push_back(curr_com[1]);
			//rfx.push_back(dBodyGetPosition(bodies[3])[0]);
			//rfy.push_back(dBodyGetPosition(bodies[3])[1]);
			if (novelty_function % 2 ==1)
				lastdown=(-1);
			else		
				lastdown=0;  //if this is set to 0, we don't care if feet sequence alternates
		}
		/*
		if(step==1 || (step-lft[lft.size()-1] > 200 && lastdown!=(-1)))
		{
			lft.push_back(step);
			lfx.push_back(dBodyGetPosition(bodies[3])[0]);
			lfy.push_back(dBodyGetPosition(bodies[3])[1]);
			lastdown=(-1);
		}
		*/
	}
	//don't let first recorded instance of both feet down set the lastdown criteria
	if(step==1)
		lastdown=0;

	leftdown=newleftdown;
	rightdown=newrightdown;
	//reset ground sensors for feetz
	onground[0]=false;
	onground[3]=false;
 }
};

static dWorldID world;
static dSpaceID space;
static dGeomID floorplane;

static vector<dBodyID> bodies;
static vector<dGeomID> geoms;
static vector<Creature*> creatures;


static dJointGroupID contactgroup;

static bool show_contacts = true;

#define CYLRADIUS    0.6
#define CYLLENGTH    2.0
#define SPHERERADIUS 0.5

    
#ifdef dDOUBLE
#define dsDrawCapsule dsDrawCapsuleD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawSphere dsDrawSphereD
#define dsDrawBox dsDrawBoxD
#define dsDrawLine dsDrawLineD
#endif



// this is called by dSpaceCollide when two objects in space are
// potentially colliding.

static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
  dBodyID b1,b2;
  dBodyID test;
  assert(o1);
  assert(o2);

  b1 = dGeomGetBody(o1);
  b2 = dGeomGetBody(o2);

  if (b1 && b2 && dAreConnected (b1,b2)) return;

  if(o1 == floorplane || o2 == floorplane)
  {
	if(o1==floorplane)
		test=b2;
	if(o2==floorplane)
		test=b1;
	//test should equal the body that is colliding with floor

	for(int x=0;x<creatures.size();x++)
	{
		int bsize=creatures[x]->bodies.size();
		for(int y=0;y<bsize;y++)
			if (test==creatures[x]->bodies[y])
				creatures[x]->onground[y]=true;
	}
  }


  const int N = 32;
  dContact contact[N];
  int n = dCollide (o1,o2,N,&(contact[0].geom),sizeof(dContact));
  if (n > 0) 
  {
    for (int i=0; i<n; i++) 
    {
      contact[i].surface.mode = 0;
      contact[i].surface.mu = dInfinity; //50.0; // was: dInfinity
      dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
      dJointAttach (c, dGeomGetBody(contact[i].geom.g1), dGeomGetBody(contact[i].geom.g2));
    }
  }
}


// start simulation - set viewpoint
#ifdef GRAPHICS
static void start()
{
  static float xyz[3] = {-2,-2,2};
  static float hpr[3] = {45.0000f,-27.5000f,0.0000f};
  dsSetViewpoint (xyz,hpr);
}


// called when a key pressed

static void command (int cmd)
{
  switch (cmd) 
  {
    case ' ':
      break;
  }
}


void drawGeom (dGeomID g, const dReal *pos, const dReal *R, int show_aabb)
{
  int i;

  if (!g) return;
  if (!pos) pos = dGeomGetPosition (g);
  if (!R) R = dGeomGetRotation (g);

  int type = dGeomGetClass (g);
  if (type == dBoxClass) {
    dVector3 sides;
    dGeomBoxGetLengths (g,sides);
    dsDrawBox (pos,R,sides);
  }
  else if (type == dSphereClass) {
    dsDrawSphere (pos,R,dGeomSphereGetRadius (g));
  }
  else if (type == dCapsuleClass) {
    dReal radius,length;
    dGeomCapsuleGetParams (g,&radius,&length);
    dsDrawCapsule (pos,R,length,radius);
  }
  else if (type == dCylinderClass) {
    dReal radius,length;
    dGeomCylinderGetParams (g,&radius,&length);
    dsDrawCylinder (pos,R,length,radius);
  }
  else if (type == dGeomTransformClass) {
    dGeomID g2 = dGeomTransformGetGeom (g);
    const dReal *pos2 = dGeomGetPosition (g2);
    const dReal *R2 = dGeomGetRotation (g2);
    dVector3 actual_pos;
    dMatrix3 actual_R;
    dMULTIPLY0_331 (actual_pos,R,pos2);
    actual_pos[0] += pos[0];
    actual_pos[1] += pos[1];
    actual_pos[2] += pos[2];
    dMULTIPLY0_333 (actual_R,R,R2);
    drawGeom (g2,actual_pos,actual_R,0);
  }

  if (show_aabb) {
    // draw the bounding box for this geom
    dReal aabb[6];
    dGeomGetAABB (g,aabb);
    dVector3 bbpos;
    for (i=0; i<3; i++) bbpos[i] = 0.5*(aabb[i*2] + aabb[i*2+1]);
    dVector3 bbsides;
    for (i=0; i<3; i++) bbsides[i] = aabb[i*2+1] - aabb[i*2];
    dMatrix3 RI;
    dRSetIdentity (RI);
    dsSetColorAlpha (1,0,0,0.5);
    dsDrawBox (bbpos,RI,bbsides);
  }
}
// render thingee

  void render_body(dBodyID body, dGeomID geom)
  {
	  const dReal *CPos = dBodyGetPosition(body);
	  const dReal *CRot = dBodyGetRotation(body);
	  dReal cpos[3] = {CPos[0], CPos[1], CPos[2]};
	  dReal crot[12] = { CRot[0], CRot[1], CRot[2], CRot[3], CRot[4], CRot[5], CRot[6], CRot[7], CRot[8], CRot[9], CRot[10], CRot[11] };

	  if(dGeomGetClass(geom)==dCylinderClass)
	  {
		  dReal rad,length;
		  dGeomCylinderGetParams(geom,&rad,&length);
   	      dsDrawCylinder
		  (
			cpos,
			crot,
			length,
			rad
		  ); // single precision
	  }
	  else if(dGeomGetClass(geom)==dSphereClass)
	  {
		  dReal rad=dGeomSphereGetRadius(geom);

		    dsDrawSphere
		  (
			cpos,
			crot,
			rad
		  ); // single precision
	  }
	  else if(dGeomGetClass(geom)==dBoxClass)
	  {
		dVector3 sides;
		dGeomBoxGetLengths(geom,sides);
		dsDrawBox
			(
			cpos,
			crot,
			sides
			);
	  }
 }
#endif
// simulation loop

void simulationStep()
{
	double timestep=0.01;
	if(!bMoviePlay)
	{
		dSpaceCollide (space,0,&nearCallback);
		dWorldStep(world,timestep);
	}
	for(int x=0;x<creatures.size();x++)
		  creatures[x]->Update(timestep);

	if(!bMoviePlay)
		dJointGroupEmpty (contactgroup);
}

void create_world(Controller* controller,bool log=false)
{
 // create world
  dRandSetSeed(10);
  dInitODE();
  
  world = dWorldCreate();
  space = dHashSpaceCreate (0);
  contactgroup = dJointGroupCreate (0);
  dWorldSetGravity (world,0,0,-9.8);
  floorplane = dCreatePlane (space,0,0,1, 0.0);
  dWorldSetERP(world,0.1);
  dWorldSetCFM(world,1E-4);
  
	Biped* biped = new Biped(log,bMoviePlay);
	dVector3 pos={0.0,0.0,0.0};

	biped->Create(world,space,pos,controller);
	//cout << "Total mass:" << biped->TotalMass() << endl;
	creatures.push_back(biped);
}

void destroy_world()
{
  dJointGroupEmpty (contactgroup);
  dJointGroupDestroy (contactgroup);

for(int x=0;x<geoms.size();x++)
  dGeomDestroy(geoms[x]);

for(int x=0;x<creatures.size();x++)
{
	creatures[x]->Destroy();
	delete creatures[x];
}
creatures.clear();
bodies.clear();
geoms.clear();

  dSpaceDestroy (space);
  dWorldDestroy (world);
  dCloseODE();
}

#ifdef GRAPHICS
static void simLoop (int pause)
{
  static int timestep=0;
	
  if(bDisplay && (creatures[0]->abort()||timestep>1000) && bDoEvolution)
	{
		cout <<"Teriminating walker..." << endl;
		((Biped*)creatures[0])->print_behavior();
		bDisplay=false;
		timestep=0;
		destroy_world();
		delete controller;
	}
	
  if(!bDisplay)
  {
	
	if(1)
	{
	  biped_epoch(neatpop,Novelty);
	  
	NEAT::Organism* champ;
	double max_fitness=0;
	std::vector<NEAT::Organism*>::iterator curorg;
	champ=*(neatpop->organisms.begin()); //Make sure at least something is chosen
	//Find the population champ
	for(curorg = neatpop->organisms.begin(); curorg != neatpop->organisms.end(); ++curorg) {
		if (((*curorg)->fitness>max_fitness)) {
			champ=(*curorg);
			max_fitness=champ->fitness;
		}
	}
	  // Generate substrate genome from input file
	  char curword[100];
	  int id;
	  ifstream substrate_file(substrate_fn, ios::in);
	  substrate_file >> curword >> id;
	  Genome* substrate_genome = new Genome(id, substrate_file);

	  // Create controller using generated substrate from champ CPPN
	  controller = new CTRNNController(substrate_genome->genesis(id, champ->net));
	  //controller = new CTRNNController(champ->net);

	  delete substrate_genome;
	}
	create_world(controller,true);
	cout << "Walker created..." << endl;
	bDisplay=true;
  }

  if (!pause) //&& !creatures[0]->abort())
  {
	 if(bSlowDown)
		usleep(10000);
	  simulationStep();
	  timestep+=1;  
  
	 if(timestep%100==0)
		cout << creatures[0]->fitness() << endl;
  } 

  dsSetColorAlpha (0.3,1.0,1.0,1.0);
  dsSetTexture(DS_WOOD);

  for(int x=0;x<geoms.size();x++)
  {
	drawGeom(geoms[x],0,0,false);
  }


  for(int x=0;x<creatures.size();x++)
  {
	  for(int y=0;y<creatures[x]->geoms.size();y++)
		drawGeom(creatures[x]->geoms[y],0,0,false);
  }
}
#endif
void update_stuff(vector<float> &k, Creature* c,bool good=true)
{

if(novelty_function==NF_FITSAMP || novelty_function==NF_FIT)
	if(good)
	k.push_back(c->fitness());
	else
	k.push_back(0.0);

if(novelty_function==NF_FITSQ || novelty_function==NF_FITSQSAMP)
{
	double f=c->fitness();
	if(good)
	k.push_back(f*f);
	else
	k.push_back(0.0);
}
if(novelty_function==NF_FITCU || novelty_function==NF_FITCUSAMP)
{
	double f=c->fitness();
	if(good)
	k.push_back(f*f*f);
	else
	k.push_back(0.0);
}

if(novelty_function==NF_COGSAMP || novelty_function == NF_COG || novelty_function==NF_COGSAMPSQ || novelty_function== NF_COGSQ || novelty_function == NF_COGCU || novelty_function==NF_COGCUSAMP)
			{
				dVector3& o_com= ((Biped*)c)->orig_com;
				dVector3& c_com= ((Biped*)c)->curr_com;
				dVector3 com;
				dVector3 delta;
				c->CenterOfMass(com);
				//calculate_delta(o_com,c_com,delta);
				calculate_delta(o_com,com,delta);
				if(novelty_function==NF_COGSAMPSQ || 
					novelty_function == NF_COGSQ)
					calculate_power(delta,2);
				
				if(novelty_function==NF_COGCUSAMP ||
					novelty_function == NF_COGCU)
					calculate_power(delta,3);

				if(good)
				{
				k.push_back(delta[0]);
				k.push_back(delta[1]);
				}
				else
				{
				k.push_back(0.0);
				k.push_back(0.0);
				}
			}
			if(novelty_function==NF_LEGSAMP || novelty_function==NF_LEG || novelty_function==NF_LEGSAMPSQ || novelty_function == NF_LEGSQ || novelty_function==NF_LEGCU || novelty_function==NF_LEGCUSAMP)
			{
				Biped* b=(Biped*)c;
				dVector3& oleft=((Biped*)c)->orig_left;
				dVector3& oright=((Biped*)c)->orig_right;
				dVector3 nleft;
				dVector3 nright;
				dVector3 ldelta;
				dVector3 rdelta;
				nleft[0]= b->lfx.back(); //dBodyGetPosition(c->bodies[0])[0];
				nleft[1]= b->lfy.back(); //dBodyGetPosition(c->bodies[0])[1];
				nleft[2]= 0.0;	
				nright[0]=b->rfx.back(); //dBodyGetPosition(c->bodies[3])[0];
				nright[1]=b->rfy.back(); //dBodyGetPosition(c->bodies[3])[1];
				nright[2]=0.0;
				calculate_delta(oleft,nleft,ldelta);
				calculate_delta(oright,nright,rdelta);
				
				if(novelty_function==NF_LEGSQ || 
					novelty_function==NF_LEGSAMPSQ)
				{
					calculate_power(ldelta,2);
					calculate_power(rdelta,2);
				}

				if(novelty_function==NF_LEGCU ||
					novelty_function==NF_LEGCUSAMP)
				{		
					calculate_power(ldelta,3);
					calculate_power(rdelta,3);
				}
				if(!good)
				{
				ldelta[0]=0.0;
				ldelta[1]=0.0;
				rdelta[0]=0.0;
				rdelta[1]=0.0;
				}
				k.push_back(ldelta[0]);
				k.push_back(ldelta[1]);
				k.push_back(rdelta[0]);
				k.push_back(rdelta[1]);		
			}
}

dReal evaluate_controller(Controller* controller,noveltyitem* ni=NULL,data_record* record=NULL,bool log=false)
{
	vector<float> k;
	dReal fitness;
	int timestep=0;
	const int simtime=1500;
	create_world(controller,log);
	while(!creatures[0]->abort() && timestep<simtime)
	{
		simulationStep();

		timestep+=1;
		
		if(timestep%100 == 0 && novelty_function % 2 == 1)
		{
			update_stuff(k,creatures[0]);
		}
		if(log && timestep%100==0)
			cout << creatures[0]->fitness() << endl;
	}
	
	if(novelty_function%2==1)
	{
		for(int x=timestep+1;x<=simtime;x++)
			if(x%100==0)
				update_stuff(k,creatures[0]); //,false);
	}
	else
	{
		update_stuff(k,creatures[0]);
	}

	//cout << timestep << endl;
	fitness=creatures[0]->fitness();
	((Biped*)creatures[0])->lft.push_back(timestep);
	((Biped*)creatures[0])->rft.push_back(timestep);
	

	//((Biped*)creatures[0])->print_behavior();

	if (ni!=NULL)
	{
		ni->novelty_scale = 1.0; //1.0 / avgDif;
				
		if(novelty_function==NF_SS)
		{
		ni->data.push_back(((Biped*)creatures[0])->lft);
		ni->data.push_back(((Biped*)creatures[0])->lfx);
		ni->data.push_back(((Biped*)creatures[0])->lfy);
		ni->data.push_back(((Biped*)creatures[0])->rft);
		ni->data.push_back(((Biped*)creatures[0])->rfx);
		ni->data.push_back(((Biped*)creatures[0])->rfy);
		}
		
		else
		{
		ni->data.push_back(k);
		}
		/*

		dVector3 com;
		creatures[0]->CenterOfMass(com);
		
		k.push_back(dBodyGetPosition(creatures[0]->bodies[0])[0]);
		k.push_back(dBodyGetPosition(creatures[0]->bodies[0])[1]);	
		k.push_back(dBodyGetPosition(creatures[0]->bodies[3])[0]);
		k.push_back(dBodyGetPosition(creatures[0]->bodies[3])[1]);
		*/
		/*
		k.push_back(com[0]);
		k.push_back(com[1]);
		*/

		//k.push_back(creatures[0]->fitness());
		
		//enable for simple novelty metrics
		//ni->data.push_back(k);
	}

	if(record!=NULL)
	{
		dVector3 com;
		creatures[0]->CenterOfMass(com);
		record->ToRec[0]=fitness;
		record->ToRec[1]=com[0];
		record->ToRec[2]=com[1];
		record->ToRec[3]=com[2];
		record->ToRec[4]=timestep;
	}

	destroy_world();
	return fitness;
}

/*
class WalkerDomain:public Domain
{
	virtual void evaluate(Individual* ind)
    {
		dReal parameters[9];
        FloatGenome* g=(FloatGenome*)ind->genome;
		
		for(int x=0;x<9;x++)
			parameters[x]=g->vals[x];

		SineController* controller = new SineController(6);
		controller->load_params(parameters);
		double fit=evaluate_controller(controller);

		ind->fitness->fitness=fit;
        ind->fitness->raw_fitness=fit;
    }
};

Individual* generator()
{
    Fitness *f=new Fitness();
    FloatGenome* g=new FloatGenome(9,-1.0,1.0);
	for(int x=0;x<9;x++)
		g->vals[x]=0.0;
	g->mutate();
    return new Individual(g,f);
}
*/

int main (int argc, char **argv)
{
  // setup pointers to drawstuff callback functions
  char curword[100];
  int id;

  #ifdef GRAPHICS
  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &start;
  fn.step = &simLoop;
  fn.command = &command;
  fn.stop = 0;
  fn.path_to_textures = DRAWSTUFF_TEXTURE_PATH;
  #endif 

  srand( (unsigned)time( NULL ) );
  load_neat_params("params.ne", true);

  //load_parameters();
  Genome* start_genome;
  Controller* controller;

  if(argc > 1)
  {
     if(strcmp(argv[1], "display") == 0)
     {
   
	    if (argc == 2)
	    {
		controller = new DummyController();
	    }
	    else
	    {
		    ifstream iFile(argv[2], ios::in);
		    ifstream substrate_file(argv[3], ios::in);
		    strcpy(substrate_fn, argv[3]);

		    cout<<"Reading in the start genome"<<endl;
		    //Read in the start Genome
		    iFile >> curword;
		    iFile >> id;
		    cout << "Reading in Genome id " << id << endl;
		    start_genome = new Genome(id, iFile);
		    iFile.close();

		    // Test activation function code RSO
		    /*cout<<"Start Genome: "<<start_genome<<endl;

		    Network* net = start_genome->genesis(0);
		    double vals[3]={0.5,0.5,1.0};
		    net->load_sensors(vals);
		    net->activate();
		    //net->activate();

		    cout << "Output 0: " << net->outputs[0]->activation << endl;
		    cin >> vals[0];*/
		    // End test activation function code RSO

		    /*
		    NEAT::Network* net =     start_genome->genesis(id);
		    net->ctrnn_dynamics();
		    delete net;
		    */
		    if (argc==3)
			    return 0;

		    // Create CPPN from start genome
		    Network* cppn_network = start_genome->genesis(id);

		    // Generate substrate genome from input file
		    cout << "Reading in the substrate genome" << endl;
		    substrate_file >> curword >> id;
		    cout << "Reading in substrate genome id " << id << endl;
		    Genome* substrate_genome = new Genome(id, substrate_file);

		    Network* substrate_network = substrate_genome->genesis(id, cppn_network);

		    for(int i = 0; i < substrate_network->all_nodes.size(); i++)
		      {
			printf("node %d, bias: %f, tc: %f", (i+1), substrate_network->all_nodes[i]->bias, substrate_network->all_nodes[i]->time_const);

			for(int j = 0; j < substrate_network->all_nodes[i]->outgoing.size(); j++)
			  {
			    printf("\n\tlink %d, from pos (%f, %f, %f), to pos (%f, %f, %f), weight: %f", (j+1),
				   substrate_network->all_nodes[i]->outgoing[j]->in_node->x_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->in_node->y_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->in_node->z_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->out_node->x_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->out_node->y_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->out_node->z_sub_pos,
				   substrate_network->all_nodes[i]->outgoing[j]->weight);
			  }
			printf("\n\n");
		      }

		    // Use CPPN to generate the substrate from the substrate genome and create
		    // the controller with the subtrate Network
		    controller = new CTRNNController(substrate_genome->genesis(id, cppn_network), true);

		    delete substrate_genome;
		    delete cppn_network;
	    }

	    bSlowDown = true;
	    bDisplay = true;
	    bDoEvolution = false;
     }
     else
	{
		strcpy(outdir, argv[2]);

		if(strcmp(argv[3], "novelty") == 0)
			Novelty = true;
		else
			Novelty = false;

		if (argc >= 5)
		{
			strcpy(startgenes_fn, argv[4]);
			cout << "using " << startgenes_fn << " for start genes" << endl;
		}

		if (argc >= 6)
		{
			strcpy(substrate_fn, argv[5]);
			cout << "using " << substrate_fn << " for substrate genes" << endl;
		}

		if (argc >= 7)
		{
			novelty_function = atoi(argv[6]);
		}

		if (argc >= 8)
		{
			srand(atol(argv[7]));
		}
	}	
  }
  

 if (bDoEvolution)
   {
     neatpop=biped_realtime(Novelty);
     if(bNoVis)
       {
	 for(int k=0; k < 20001; k++)
	   {
	     biped_epoch(neatpop, Novelty);
	}    
       }
   }

//  neatpop->organisms[0]->net->ctrnn_dynamics();
//  return 0;

/*
if(0)
	{
	domain = new WalkerDomain();
	pop = new Population();
	ga= new SteadyGA(pop,domain,&generator,100);
	}
	else
	{
		SineController* controller = new SineController(6);
		controller->load_params(params);
		create_world(controller,true);
		bDisplay=true;
		bDoEvolution=false;
	}
*/

// run simulation
#ifdef GRAPHICS
 if(!bNoVis)
   {
     if(!bDoEvolution)
       create_world(controller,true);
     dsSimulationLoop (argc,argv,352,288,&fn);
     destroy_world();
   }
#endif
 cout << bNoVis << endl;
 if(bNoVis && !bDoEvolution)
   {
     cout << "EVALUATING CONTROLLER..." << endl;
     cout << evaluate_controller(controller,NULL,NULL,true) << endl;
     ((Biped*)creatures[0])->print_behavior();
   }
  //delete neatpop;
  return 0;
}

NEAT::Population *biped_realtime(bool novelty) {
	NEAT::Population *pop;
	NEAT::Genome *start_genome;
    char curword[20];
    int id;

    ostringstream *fnamebuf;
    int gen;

    double highscore;

    ifstream iFile(startgenes_fn,ios::in);

    cout<<"Reading in the start genome"<<endl;
    //Read in the start Genome
    iFile>>curword;
    iFile>>id;
    cout<<"Reading in Genome id "<<id<<endl;
	start_genome=new NEAT::Genome(id,iFile);
    iFile.close();

    cout<<"Start Genome: "<<start_genome<<endl;

    //Spawn the Population from starter gene
    cout<<"Spawning Population off Genome"<<endl;
	pop=new NEAT::Population(start_genome,NEAT::pop_size);
      
    //Alternative way to start off of randomly connected genomes
    //pop=new Population(pop_size,7,1,10,false,0.3);

    cout<<"Verifying Spawned Pop"<<endl;
    pop->verify();
    
    //Start the evolution loop using rtNEAT method calls 
    biped_realtime_loop(pop, novelty);
    
    delete start_genome;

    return pop;
}

void biped_epoch(NEAT::Population *pop,bool novelty) {

vector<NEAT::Organism*>::iterator curorg;
	vector<NEAT::Species*>::iterator curspecies;

	vector<NEAT::Species*>::iterator curspec; //used in printing out debug info                                                         

	vector<NEAT::Species*> sorted_species;  //Species sorted by max fit org in Species                                                  

  static int epoch=0;
  int pause;
  bool win=false;

  double champ_fitness;
  NEAT::Organism *champ;

  //Real-time evolution variables                                                                                             
  int offspring_count;
  NEAT::Organism *new_org;

  //We try to keep the number of species constant at this number                                                    
  int num_species_target=NEAT::pop_size/15;
  
  //This is where we determine the frequency of compatibility threshold adjustment
  int compat_adjust_frequency = NEAT::pop_size/10;
  if (compat_adjust_frequency < 1)
    compat_adjust_frequency = 1;


  //Rank all the organisms from best to worst in each species
  pop->rank_within_species();                                                                            

  //Assign each species an average fitness 
  //This average must be kept up-to-date by rtNEAT in order to select species probabailistically for reproduction
  pop->estimate_all_averages();
  
  epoch++;

  //if(epoch>1) exit(1);
  if(epoch%2500==0)
  {
	char file[50];
	
	//sprintf(file,"%sgen%d",outdir,epoch/250);
	//pop->print_to_file_by_species(file);
	
	//sprintf(file,"%sarchive.dat",outdir);
	//archive.Serialize(file);

	sprintf(file,"%srecord.dat",outdir);
	Record.serialize(file);

        sprintf(file,"%sfittest",outdir);
	archive.serialize_fittest(file);
  }
  


  //Now create offspring one at a time, testing each offspring,                                                               
  // and replacing the worst with the new offspring if its better
  for (offspring_count=0;offspring_count<250;offspring_count++) {
        
    //Every pop_size reproductions, adjust the compat_thresh to better match the num_species_targer
    //and reassign the population to new species                                              
    if (offspring_count % compat_adjust_frequency == 0) {

	if(novelty)
	{	
	   //update fittest individual list		
	   archive.update_fittest(pop);
	   //refresh generation's novelty scores
	   archive.evaluate_population(pop,true);
	}
      int num_species = pop->species.size();
      double compat_mod=0.1;  //Modify compat thresh to control speciation                                                     

      // This tinkers with the compatibility threshold 
      if (num_species < num_species_target) {
	NEAT::compat_threshold -= compat_mod;
      }
      else if (num_species > num_species_target)
	NEAT::compat_threshold += compat_mod;

      if (NEAT::compat_threshold < 0.3)
	NEAT::compat_threshold = 0.3;

      cout<<"compat_thresh = "<<NEAT::compat_threshold<<endl;

      //Go through entire population, reassigning organisms to new species                                                  
      for (curorg = (pop->organisms).begin(); curorg != pop->organisms.end(); ++curorg) {
	pop->reassign_species(*curorg);
      }
    }
    
/*
    //For printing only
    for(curspec=(pop->species).begin();curspec!=(pop->species).end();curspec++) {
      cout<<"Species "<<(*curspec)->id<<" size"<<(*curspec)->organisms.size()<<" average= "<<(*curspec)->average_est<<endl;
    }

    cout<<"Pop size: "<<pop->organisms.size()<<endl;
*/

    //Here we call two rtNEAT calls: 
    //1) choose_parent_species() decides which species should produce the next offspring
    //2) reproduce_one(...) creates a single offspring fromt the chosen species
    new_org=(pop->choose_parent_species())->reproduce_one(offspring_count,pop,pop->species);

    //Now we evaluate the new individual
    //Note that in a true real-time simulation, evaluation would be happening to all individuals at all times.
    //That is, this call would not appear here in a true online simulation.
    //cout<<"Evaluating new baby: "<<endl;

	data_record* newrec=new data_record();
	newrec->indiv_number=indiv_counter;
		
	new_org->noveltypoint = biped_evaluate(new_org,newrec);
	new_org->noveltypoint->indiv_number = indiv_counter;
	//calculate novelty of new individual
	if(novelty)
	{
	archive.evaluate_individual(new_org,pop);
	new_org->fitness*=new_org->noveltypoint->novelty_scale;
	archive.update_fittest(new_org);
	}	
	else
	{
	new_org->fitness=new_org->noveltypoint->fitness;
	}
	
	//add record of new indivdual to storage
	indiv_counter++;
	Record.add_new(newrec);
	
	//update fittest list

    if (win) {
      cout<<"WINNER"<<endl;
      //pop->print_to_file_by_species("rt_winpop");
      break;
    }

    //Now we reestimate the baby's species' fitness
    new_org->species->estimate_average();

    //Remove the worst organism                                                                                               
    pop->remove_worst();

  }
	
  if(novelty)
    {
      archive.end_of_gen_steady(pop);
      //archive.add_randomly(pop);
      archive.evaluate_population(pop,false);
      cout << "ARCHIVE SIZE:" << 
	archive.get_set_size() << endl;
    }
}


void biped_realtime_loop(NEAT::Population *pop,bool novelty) {
  vector<NEAT::Organism*>::iterator curorg;
  vector<NEAT::Species*>::iterator curspecies;
  
  vector<NEAT::Species*>::iterator curspec; //used in printing out debug info                                                         
  
  vector<NEAT::Species*> sorted_species;  //Species sorted by max fit org in Species                                                  

  int pause;
  bool win=false;

  double champ_fitness;
  NEAT::Organism *champ;

  //Real-time evolution variables                                                                                             
  int offspring_count;
  NEAT::Organism *new_org;

  //We try to keep the number of species constant at this number                                                    
  int num_species_target=NEAT::pop_size/15;
  
  //This is where we determine the frequency of compatibility threshold adjustment
  int compat_adjust_frequency = NEAT::pop_size/10;
  if (compat_adjust_frequency < 1)
    compat_adjust_frequency = 1;
  
  //Initially, we evaluate the whole population                                                                               
  //Evaluate each organism on a test                                                                                          
  for(curorg=(pop->organisms).begin();curorg!=(pop->organisms).end();++curorg)
    {
      
      //shouldn't happen                                                                                                        
      if (((*curorg)->gnome)==0) {
	cout<<"ERROR EMPTY GENOME!"<<endl;
	cin>>pause;
      }
      
      (*curorg)->noveltypoint=(biped_evaluate((*curorg)));
      (*curorg)->noveltypoint->indiv_number=indiv_counter;
      (*curorg)->fitness=(*curorg)->noveltypoint->fitness;
      indiv_counter++;
    }

  //Get ready for real-time loop
  if(novelty)
    {
      //assign fitness scores based on novelty
      archive.evaluate_population(pop,true);
      //add to archive
      archive.evaluate_population(pop,false);
    }

}

noveltyitem* biped_evaluate(NEAT::Organism *org,data_record* data)
{
	char curword[100];
	int id = 0;

	noveltyitem *new_item = new noveltyitem;
	
	// Generate substrate genome from input file
	ifstream substrate_file(substrate_fn, ios::in);
	substrate_file >> curword >> id;
	
	Genome* substrate_genome = new Genome(id, substrate_file);
	
	// Store substrate's genome as the new novelty item's genotype
	//new_item->genotype = substrate_genome;
	new_item->genotype = new Genome(*org->gnome);

	// Create CPPN from the Organism's genome
	Network* cppn_network = org->gnome->genesis((*org->net).net_id);
	//Network* cppn_network = org->net;
	
	// Use CPPN to generate the substrate from the substrate genome and store
	// the substrate network as the phenotype
	Network* substrate_network = substrate_genome->genesis((*org->net).net_id, cppn_network);
  	new_item->phenotype = cppn_network;

	// test code - prints out generated substrate
	/*for(int i = 0; i < substrate_network->all_nodes.size(); i++)
	{
		printf("node %d, bias: %f, tc: %f", (i+1), substrate_network->all_nodes[i]->bias, substrate_network->all_nodes[i]->time_const);

		for(int j = 0; j < substrate_network->all_nodes[i]->outgoing.size(); j++)
		{
			printf("\n\tlink %d, from pos (%f, %f, %f), to pos (%f, %f, %f), weight: %f", (j+1),
			substrate_network->all_nodes[i]->outgoing[j]->in_node->x_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->in_node->y_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->in_node->z_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->out_node->x_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->out_node->y_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->out_node->z_sub_pos,
			substrate_network->all_nodes[i]->outgoing[j]->weight);
		}
		printf("\n\n");
	}
	scanf("%d", &id);*/
	// end test code
	
	// Create the controller with the substrate Network
	CTRNNController* cont = new CTRNNController(substrate_network);
	new_item->fitness = evaluate_controller(cont, new_item, data);

	// Clean up objects that aren't needed any more
	delete cont;
	delete substrate_network;
	delete substrate_genome;

	return new_item;
}

// Scales given val in between the given min and max
double NEAT::scale(double val, double min, double max)
{
	double scaled = val;

	if(scaled < -1.0)
		scaled = -1.0;

	else if(scaled > 1.0)
		scaled = 1.0;

	scaled = ((scaled + 1.0) / 2.0) * (max - min) + min;

	return scaled;
}
