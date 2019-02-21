
var mass_spring_system;
var flock;

var canvas_size = 600;
const rad_fac = 0.70; // sets radius of bubble
const dx = 1;  // grid spacing
const dt = 0.3; // timestep

const nnodes  = 80; // numbers of mass nodes
const node_mass = 6; // mass of nodes
const nboids = 140;
const boid_mass = 0.3;

var ks = 3/dx; // spring constant
const gammas = 0.0;   // spring damping parm

const gamma_node = 0.004;  // damping on nodes (force depends on velocity)

const force_amp = 10.0; // for interactions between boids and nodes
const force_k = 0.2;  // 1/scale for interactions between boids and nodes
const vforce_amp = 0.00;  // damping


// distances for boid forces 
const d_repel = 50;
const d_attract = 20;  // scale over which attractive force applies
const d_align = 50; // scale for alignment
// strengths for boid forces (accelerations)
const attract_force = 0.00*boid_mass; // amplitude of attract force
const repel_force = 0.08*boid_mass;  // repel
const align_force = 0.05*boid_mass;

const boidspeed = 4;  // speed of self propulsion

const ndt=3;  // number of time steps per display update

const propel_force = 1; // not currently used

function setup() {
  createCanvas(canvas_size, canvas_size); 
  let xwidth = width*dx;     // grid size set here!
  let yheight = height*dx;
  let big_radius = rad_fac*xwidth/2;
  // set up masses/springs
  mass_spring_system = new Mass_spring_system(nnodes,big_radius);
    
  mass_spring_system.display_springs();
  mass_spring_system.display_nodes();
  
  // set up flock  of boids
  flock = new Flock(nboids);
  flock.display_flock();
  // let boid_set = flock.boid_set;
  // let node_set = mass_spring_system.node_set;
  // boid_node_interact(boid_set,node_set,force_amp,force_k,vforce_amp);
  
}

var dcount=0; // used for centroiding display

function draw() {

   let boid_set = flock.boid_set;
   let node_set = mass_spring_system.node_set;
   background(240);
    
   mass_spring_system.display_springs();
   mass_spring_system.display_nodes();
   flock.display_flock();
    
   for (let k=0;k<ndt;k++){ // numbers of timesteps per display
     
     // zero accelerations
      mass_spring_system.zeroaccel();
      flock.zeroaccel();   
    
      // compute accelerations on mass/spring system
      mass_spring_system.compute_accel();
      // compute accelerations on boid flock
      // flock.propel();
      flock.align();
      flock.cohesion();
      flock.repel();
      // compute interactions boids/masses
      boid_node_interact(boid_set,node_set,force_amp,force_k,vforce_amp);
    
      // update
      mass_spring_system.single_timestep(dt);
      flock.single_timestep(dt);
  
      if (nnodes<2){
        flock.borders(); // if there aren't any mass nodes
      }
      else{ // centroid display
        dcount++;
        if ((dcount%10)==0){  // shift centroid
          let centroid = mass_spring_system.centroid();
          mass_spring_system.shift(centroid);
          flock.shift(centroid);
        }
      }
   }
}

function Flock(nboids) {
  // An array for all the boids
  this.boid_set = []; // Initialize the array of boids
  this.dx = dx;
  // this.maxforce = 20;
  // this.maxspeed = 2;
    
  this.boidspeed = boidspeed;
  
  this.d_attract   = d_attract; // distance for attraction/repulsion forces between Boids
  this.d_repel   = d_repel;  // distance for attraction/repulsion forces between Boid
  this.d_align   = d_align;  // distance for alignment
  
  this.attract_force = attract_force;  // attact force size
  this.repel_force   = repel_force;    // repel force size
  this.align_force   = align_force;
  this.propel_force  = propel_force;

  
  // create nboids of Boids!
  for (let i=0;i<nboids;i++){
    let b = new Boid(this.dx);
    this.boid_set.push(b);
  }
  
  this.display_flock = function(){  // display boids
    let n = this.boid_set.length;
    for (let i=0;i<n;i++){
      let b = this.boid_set[i];
      b.render();
    }
  }
  this.borders = function(){
    let n = this.boid_set.length;
    for (let i=0;i<n;i++){
      let b = this.boid_set[i];
      b.borders();
    }
  }
  
  // zero all the accelerations
  this.zeroaccel = function(){
        let n = this.boid_set.length;
    for(let k = 0;k<n;k++){
      this.boid_set[k].acceleration.mult(0);
    }
  }
  
  // Repel force depends on inverse distance -- all boid pairs
  this.repel = function() {
      let n = this.boid_set.length;
      for (let i = 0; i < n-1; i++) {
      let bi = this.boid_set[i];
      for (let j = i+1; j <n; j++) {
        let bj = this.boid_set[j];
        let dr = p5.Vector.sub(bi.position,bj.position);
        let r_len = dr.mag();  // length of interboid distance
        if (r_len < this.d_repel){
          let drhat = dr.copy();
          drhat = drhat.normalize();  // unit vector for direction
          let Force = drhat.copy();
          let fac = -1.0*this.repel_force*this.d_repel/r_len; // normalized here
          Force.mult(fac); 
          let ai = p5.Vector.mult(Force,-1/bi.m);
          let aj = p5.Vector.mult(Force, 1/bj.m);
          bi.acceleration.add(ai);
          bj.acceleration.add(aj);
        }
      // down side of this method is if we have many nearby particles, acceleration gets high
      }
    }
  }
    
  // Cohesion force, doing all pairs
  this.cohesion = function() {
     if (this.attract_force >0){
     let n = this.boid_set.length;
     for (let i = 0; i < n; i++) {
        let bi = this.boid_set[i];
        let pos_sum = createVector(0,0);
        let count = 0;
        for (let j = 0; j < n; j++) {
           let bj = this.boid_set[j];
           let dr = p5.Vector.sub(bi.position,bj.position);
           let r_len = dr.mag();  // length of interboid distance
           if ((r_len > 0) && (r_len < this.d_attract)){
              pos_sum.add(bj.position); // sum of positions
              count++;
              }
           }
           if (count >0){
              pos_sum.div(count);   // is average of nearby boid positions!
              let Force = p5.Vector.sub(pos_sum, bi.position); // target direction
              Force.normalize();
              Force.mult(this.boidspeed);
              Force.sub(bi.velocity);
              let ai = p5.Vector.mult(Force,this.attract_force/bi.m);  
              bi.acceleration.add(ai);
           }
        }
     }
  }
  
  // alignment try to stear toward mean velocity of nearby Boids
  // velocity dependent forces here!
  this.align = function() {  
      let n = this.boid_set.length;
      // For every boid in the system, check if it's close to another
      for (let i = 0; i < n; i++) {
      let count = 0;
      let bi = this.boid_set[i];
      // let steer = createVector(0,0);  // from average of nearest neighbor velocities
      let v_ave = createVector(0,0);
      for (let j = 0; j < n; j++) {
        let bj = this.boid_set[j];
        let dr = p5.Vector.dist(bi.position,bj.position);
        // If the distance is greater than 0 and less than an arbitrary amount (0 when you are yourself)
        if ((dr > 0) && (dr < this.d_align)) {
            v_ave.add(bj.velocity);  // sum of velocities of nearby
            count++; // Keep track of how many nearby
        }
      }
      if (count > 0){
        v_ave.normalize();
        v_ave.mult(this.boidspeed);
        let steer = p5.Vector.sub(v_ave,bi.velocity);
        // desired is now the stear 
        // As long as the vector is greater than 0
        if (steer.mag() > 0) {
           steer.mult(this.align_force/bi.m);
           // steer.limit(this.maxforce);
           bi.acceleration.add(steer);//  only on bi
         }
      }
    }
  }
  
  // try to reach same velocity, for is propto propel_force times 
  // difference in velocity from boidspeed
  this.propel = function(){
    let n = this.boid_set.length;
    for (let i = 0; i < n; i++) {
      let bi = this.boid_set[i];
      let vi = bi.velocity.copy();
      let vmag = vi.mag();  // length
      let vhat = vi.copy(); // unit vector
      vhat.normalize();
      let Force = vhat.copy(); // direction same as velocity
      Force.mult((vmag-this.boidspeed)*this.propel_force);
      let ai = p5.Vector.div(Force,-bi.m);
      bi.acceleration.add(ai);
    }
  }
  
  // integrate dt
  this.single_timestep = function(dt){
    // this.zeroaccel();
    // this.propel(); // propel boids
    // this.repel(); // repel repulsion
    // this.cohesion(); // attraction repulsion
    // this.align();
    // boid_node_interact(boid_set,node_set,
        //        force_amp,force_k,vforce_amp); 
    // apply interactions between boids and nodes
    let n = this.boid_set.length;
    for (let i = 0; i < n; i++) {
      bi = this.boid_set[i];
      let dv = p5.Vector.mult(bi.acceleration,dt); // h
      bi.velocity.add(dv); // 
      let dr = p5.Vector.mult(bi.velocity,dt);
      bi.position.add(dr); //
    }
  }
  
    
  this.shift = function(centroid){
        let n = this.boid_set.length;
        for(let i = 0;i< n;i++){
            boidi = this.boid_set[i];
            boidi.position.sub(centroid);
        }
  }
  
}


// exponential short range forces between boids and nodes
// U = force_amp*exp(-force_k*d) where d = distance between
function boid_node_interact(boid_set,node_set,
                             force_amp,force_k,vforce_amp){
  let n_nod = node_set.length;
  let n_bd = boid_set.length;
  
  for (let i=0;i<n_nod;i++){// loop over nodes
    let nodei = node_set[i];
    for(let j=0;j<n_bd;j++){ // loop over boids
      boidj = boid_set[j];
      let dr = p5.Vector.sub(boidj.position,nodei.position);  // vector between
      let d = dr.mag(); // distance between 
      let dv = p5.Vector.sub(boidj.velocity,nodei.velocity);
      //let v = dv.mag;
      if (d < 3/force_k){
        let drhat = dr.copy();
        drhat.normalize(); //  is a unit vector
        let Force = drhat.copy(); // force in direction between
        Force.mult(-force_amp*exp(-force_k*d));
        let vForce = dv.copy();
        vForce.mult(vforce_amp);// damping force depends on velocity diff
        Force.add(vForce);
        let a_boid = Force.copy();
        let a_node = Force.copy();
        a_boid.div(-boidj.m); // acceleration
        a_node.div( nodei.m); // acceleration
        boidj.acceleration.add(a_boid);
        nodei.acceleration.add(a_node);
        
      }
    }
  }
  
}
