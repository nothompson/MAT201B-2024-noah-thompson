#include "al/app/al_App.hpp"
#include "al/app/al_GUIDomain.hpp"
#include "al/math/al_Random.hpp"
#include "Gamma/SamplePlayer.h"
#include "Gamma/Effects.h"
#include "Gamma/Analysis.h"
#include "Gamma/Envelope.h"

using namespace al;
using namespace std;

#include <fstream>
#include <vector>

#include "al/app/al_DistributedApp.hpp"
#include "al_ext/statedistribution/al_CuttleboneDomain.hpp"
#include "al_ext/statedistribution/al_CuttleboneStateSimulationDomain.hpp"


struct CommonState {
  // float particlePositions[10000];
  float colorState[3];

};



Vec3f randomVec3f(float scale) {
  return Vec3f(rnd::uniformS(), rnd::uniformS(), rnd::uniformS()) * scale;
}
string slurp(string fileName); 

struct FilterFollow {
  gam::Biquad<> filter;
  gam::EnvFollow<> follow;
  float operator()(float in) {
    return gam::scl::mapLin(gam::scl::ampTodB(follow(filter(in))), -134.0f, 0.0f, 0.0f, 1.0f);
  }
};

struct MyApp : DistributedAppWithState<CommonState> {

  Parameter parameter[3]{
    {"low", "", 0,  0, 127},
    {"mid", "", 0,  0, 127},
    {"high", "", 0, 0, 127}
  };
  FilterFollow follow[3];

  colorState[0] = parameter[0];
  colorState[1] = parameter[1];
  colorState[2] = parameter[2];

  Parameter value {"/value", "", 0,0,1};
  gam::EnvFollow<>ampFollow;
  gam::Accum<>tmr;
  gam::AD<>env;
  float slope;




  Parameter AUDIO{"/AUDIO", "", 0, 0, 0};
  Parameter amplitude{"/amplitude", "", 0.8, 0.0, 1.0};
  Parameter impulseRate{"/impulseRate", "", 1.0, 0.1, 100.0};

  Parameter rateToggle{"/rateToggle", "", 0.0, 0.0, 1.0};
  Parameter grainRateMax{"/grainRateMax", "", 1.0, 0.1, 20.0};
  Parameter grainRateMin{"/grainRateMin", "", 1.0, 0.01, 20.0};

  Parameter grainPosition{"/grainPosition", "", 0.0, 0.0, 1.0};
  Parameter grainSpread{"/grainSpread", "", 0.0, 0.0, 1.0};
  // Parameter grainPositionMin{"/grainPositionMin", "", 1.0, 0.0, 1.0};
  // Parameter grainPositionMax{"/grainPositionMax", "", 1.0, 0.0, 1.0};

  Parameter reverseAmount{"/reverseAmount", "", 0.0, 0.0, 1.0};
  // Parameter fadeIn{"/fadeIn", "", 2.0, 2.0, 10000.0};
  // Parameter fadeOut{"/fadeOut", "", 2.0, 2.0, 10000.0};

  // Parameter fadeRate{"/fadeRate", "", 1.2, 0.1, 10.0};
  Parameter fadeSlope{"/fadeSlope", "", 0.3, 0.0, 1.0};

  Parameter VISUALS{"/VISUALS", "", 0, 0, 0};
  // Parameter pointSize{"/pointSize", "", 1.0, 0.0, 2.0};
  // Parameter timeStep{"/timeStep", "", 0.1, 0.01, 0.6};
  Parameter dragFactor{"/dragFactor", "", 0.3, 0.0, 0.9};
  Parameter sphereRadius{"/sphereRadius", "", 0.5, 0.05, 1.0};
  Parameter hookeConstant{"/hookeConstant", "", 0.1, 0.0, 1.0};
  // Parameter coulombConstant{"/coulombConstant", "", 0.0, 0.0, 0.5};
  // Parameter symmetry{"asymmetry", "", 1.0, 0.01, 1.0};
  // Parameter kickAmount{"/kickAmount", "", 1.0, 0.05, 10.0};

  ShaderProgram pointShader;

  //  simulation state
  Mesh mesh;  // position *is inside the mesh* mesh.vertices() are the positions
  vector<Vec3f> velocity;
  vector<Vec3f> force;
  // vector<Vec3f> displacement;
  vector<float> mass;
  
  // Parameter panRange{"/panRange", "", 0.5, 0.0, 1.0};

  // gam::Pan<>mPan;

  

  gam::SamplePlayer<float, gam::ipl::Cubic, gam::phsInc::Loop> player;

  void onInit() override {
    player.load("growth012.wav");

    follow[0].filter.type(gam::LOW_PASS);
    follow[1].filter.type(gam::BAND_PASS);
    follow[2].filter.type(gam::HIGH_PASS);
    follow[0].filter.freq(200);
    follow[1].filter.freq(2000);
    follow[2].filter.freq(8000);

    auto GUIdomain = GUIDomain::enableGUI(defaultWindowDomain());
    auto &gui = GUIdomain->newGUI();
    gui.add(AUDIO);
    gui.add(amplitude);
    gui.add(impulseRate);
    gui.add(rateToggle);
    gui.add(grainRateMax);
    gui.add(grainRateMin);
    gui.add(grainPosition);
    gui.add(grainSpread);
    // gui.add(grainPositionMin);
    // gui.add(grainPositionMax);
    gui.add(reverseAmount);
    // gui.add(fadeRate);
    gui.add(fadeSlope);
    // gui.add(fadeIn);
    // gui.add(fadeOut);
    gui.add(VISUALS);
    // gui.add(pointSize);  // add parameter to GUI
    // gui.add(timeStep);   // add parameter to GUI
    gui.add(dragFactor);   // add parameter to GUI
    gui.add(sphereRadius);
    gui.add(hookeConstant);
    // gui.add(coulombConstant);
    // gui.add(symmetry);
    // gui.add(kickAmount);
    // gui.add(panRange);
  }

  // MyApp(){
  //   tmr.period(1 / impulseRate);
  //   tmr.phaseMax();
  //   slope = fadeSlope;
  // }

  float pointSize = (grainRateMax + grainRateMin) / 2;

  void onCreate() override {
    // compile shaders
    pointShader.compile(slurp("../point-vertex.glsl"),
                        slurp("../point-fragment.glsl"),
                        slurp("../point-geometry.glsl"));

    auto randomColor = []() { return HSV(rnd::uniform(), 1.0f, 1.0f); };

    mesh.primitive(Mesh::POINTS);
    for (int _ = 0; _ < 10000; _++) {
      mesh.vertex(Vec3f(2 * sin(1.5), 2 * cos(1.5), (1/ (2 * M_PI)) * 1.5));
      // mesh.vertex(state().particlePositions[_]);
      mesh.color(randomColor());

      // float m = rnd::uniform(3.0, 0.5);
      float m = 3 + rnd::normal() / 2;
      if (m < 0.5) m = 0.5;
      mass.push_back(m);

      // using a simplified volume/size relationship
      mesh.texCoord(pow(m, 1.0f / 3), 0);  // s, t

      // separate state arrays
      velocity.push_back(randomVec3f(0.1));
      force.push_back(randomVec3f(1));
    }

    nav().pos(0, 0, 20);

    ampFollow.lag(0.1);
  }

  bool freeze = false;
  float time = 0;
  float delay = 0;
  void onAnimate(double dt) override {
    // particle code
    if (freeze) return;

    // Octree tree(Vec3f(0), Vec3f(1), 0.25);
    // tree.build(mesh.vertices());

    float mSphere = sphereRadius + value;

    for (int i = 0; i < velocity.size(); i++) {
      Vec3f X0 = mesh.vertices()[i]; //resting 
      Vec3f X = (X0 / X0.mag()) * mSphere; //distance
      Vec3f F = -1 * hookeConstant * (X0 - X); //hooke equation
      force[i] += mass[i] * F; //implement 
    }

    // for (int i = 0; i < velocity.size(); i++) {
    //   // vector<int> indices;
    //   // tree.queryRegion(mesh.vertices()[i],Vec3f(1), indices);
    //   // for (int j= i+1; j < velocity.size();j++){
    //     for (int j: indices){
    //     // if (i !=j){ //don't check itself
    //     // Vec3f Q1 = mesh.vertices()[j]; // charge1
    //     // Vec3f Q2 = mesh.vertices()[i]; //charge 2
    //     Vec3f R = mesh.vertices()[j] - mesh.vertices()[i]; //distance 
    //     float F2 = ((mass[i] * mass[j]) / (R.mag()*R.mag())) * coulombConstant; //coulomb equation
    //     R.normalize();
    //     force[j] += (R * (mass[j] * F2) * 0.0001) * symmetry; // implement
    //     force[i] -= R * (mass[i] * F2) * 0.0001;
    //     // }
    //   }
    // }

    for (int i = 0; i < velocity.size(); i++) {
      force[i] += - velocity[i] * dragFactor;
    }

    vector<Vec3f> &position(mesh.vertices());

    float mImpulse = impulseRate;

    if (mImpulse > 10.0) {
      mImpulse = 10.0;
    }

    for (int i = 0; i < velocity.size(); i++) {
      // "semi-implicit" Euler integration
      velocity[i] += (force[i] / mass[i] * (mImpulse * 0.125)) * (value + (parameter[0]));
      position[i] += (velocity[i] * (mImpulse * 0.125));
    }


    for (auto &a : force) a.set(0);

    // granulation code

    if (time > delay) {
      time -= delay;
      // delay = (rnd::uniform(1.0, 0.5)) / impulseRate;
      delay = 1 / impulseRate;

        // for melodic samples
      if (rateToggle < 0.5){
        player.rate(grainRateMax * (rnd::prob(reverseAmount) ? -1 : 1));
      }
        // for textural samples
      if (rateToggle > 0.5){
        player.rate((float)rnd::uniform(1 * grainRateMax, 1 * grainRateMin) * (rnd::prob(reverseAmount) ? -1 : 1));
      }

      float positionMin = grainPosition - grainSpread;
      float positionMax = grainPosition + grainSpread;
      if (positionMin < 0){
        positionMin = 0;
      } 
      if (positionMax > 1){
        positionMax = 1;
      }

      player.phase((float)rnd::uniform(positionMin, positionMax));
      // player.phase((float)grainPosition);
      // player.phase(0);

    }
    time += dt;
  }

  void onSound(AudioIOData& io) override {
    // mPan.pos(rnd::uniform(1 * panRange, 0));
    tmr.phaseMax();
    slope = fadeSlope;
    while (io()) {
      // player.fade(fadeIn,fadeOut);

    tmr.period(1 / impulseRate);

      if(tmr()){
        env.attack(slope);
        env.decay(1-slope);
        env.amp(1.0);
        env.reset();

        slope +=0.01;
        if(slope > 1) slope = fadeSlope;
      }

      // float s1 = player() * env() * amplitude;
      float s1 = player() * amplitude;
      float s2 = player();
      // float s2;
      // mPan(s1, s1, s2);
      io.out(0) = s1;
      io.out(1) = s1;
      value.set(25 * ampFollow(s2));
      for (int i = 0; i < 3; i++) {
        parameter[i].set(follow[i](s2));
      }
    }

    // cout << (parameter[0]) << endl;

      // x y z forces based on 3 band filter analysis 
    for (int i = 0; i < velocity.size(); i++) {
      // force[i] += Vec3f(rnd::uniform(10,1) * parameter[0] * value * 2, rnd::uniform(10,1) * parameter[1] * value * 2, rnd::uniform(10,1)* parameter[2] * value * 2);
      // force[i] -= Vec3f(rnd::uniform(10,1) * parameter[0] * value * 2, rnd::uniform(10,1) * parameter[1] * value * 2, rnd::uniform(10,1)* parameter[2] * value * 2);
      
      force[i] += Vec3f(sin(rnd::uniform(10,1) * parameter[0] * value * 2), cos(rnd::uniform(10,1) * parameter[1] * value * 2), (rnd::uniform(10,1)* parameter[2] * value * 2) / (2 * M_PI));
      force[i] -= Vec3f(sin(rnd::uniform(10,1) * parameter[0] * value * 2), cos(rnd::uniform(10,1) * parameter[1] * value * 2), (rnd::uniform(10,1)* parameter[2] * value * 2) / (2 * M_PI));

      // force[i] += Vec3f(parameter[0] * value * 20, parameter[1] * value * 20, parameter[2] * value * 20);
      // force[i] -= Vec3f(parameter[0] * value * 20, parameter[1] * value * 20, parameter[2] * value * 20);
      
      // force[i] -= Vec3f(rnd::uniform(10,1) * parameter[0], rnd::uniform(10,1) * parameter[1], rnd::uniform(10,1)* parameter[2]);
      // force[i] += randomVec3f(value) * 5;
      
      // force[i] += randomVec3f(2 * parameter[0]);
      // force[i] += randomVec3f(2 * parameter[1]);
      // force[i] += randomVec3f(2 * parameter[2]);
    }
  }

  bool onKeyDown(const Keyboard &k) override {
    if (k.key() == ' ') {
      freeze = !freeze;
    }

    if (k.key() == '1') {
      player.load("growth012.wav");
    }

    if (k.key() == '2') {
      player.load("piano.wav");
    }

    if (k.key() == '3') {
      player.load("3patch.wav");
    }

    if (k.key() == '4') {
      player.load("myGuitarCmin..wav");
    }

    if (k.key() == '5') {
      player.load("take4.wav");
    }

  



    // if (k.key() == '1') {
    //   // introduce some "random" forces
    //   for (int i = 0; i < velocity.size(); i++) {
    //     // F = ma
    //     force[i] += randomVec3f(kickAmount);
    //   }
    // }

    return true;
  }

  void onDraw(Graphics &g) override {
    g.clear(0);
    g.shader(pointShader);
    g.shader().uniform("pointSize", pointSize / 100);
    g.blending(true);
    g.blendTrans();
    g.depthTesting(true);
    // g.rotate(360*time, 0, 1, 0);
    // g.color(RGB(0.25 + sqrt(parameter[0]), 0.25 + sqrt(parameter[1]), 0.25 + sqrt(parameter[2])));
    g.color(HSV((state().colorState[0] + state().colorState[1], state().colorState[2]),0.25 + value,1));
    g.draw(mesh);
  }
};

int main() {
  MyApp app;
  app.start();
}

string slurp(string fileName) {
  fstream file(fileName);
  string returnValue = "";
  while (file.good()) {
    string line;
    getline(file, line);
    returnValue += line + "\n";
  }
  return returnValue;
}

