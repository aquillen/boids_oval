
var mass_spring_system;
var flock;

var canvas_size = 600;
const rad_fac = 0.70; // sets radius of node boundary 
const dx = 1;  // grid spacing
const dt = 0.3; // timestep

const nnodes  = 80; // numbers of mass nodes
const node_mass = 6; // mass of nodes
const nboids = 140;  // number of boids
const boid_mass = 0.3; // boid mass

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

// sliders 
var repel_force_Tslider; 
var d_repel_Tslider;
var align_force_Tslider;
var d_align_Tslider;
var gamma_node_Tslider;
var ks_Tslider;


function setup() {
  createCanvas(canvas_size, canvas_size); 
  let xwidth = width*dx;     // grid size set here!
  let yheight = height*dx;
  let big_radius = rad_fac*xwidth/2;
  // set up masses/springs
  mass_spring_system = new Mass_spring_system(nnodes,big_radius,node_mass);
    
  mass_spring_system.display_springs();  // display springs
  mass_spring_system.display_nodes();   // display nodes
  
  // set up flock  of boids
  flock = new Flock(nboids);  
  flock.display_flock(); // display flock of boids
  // let boid_set = flock.boid_set;
  // let node_set = mass_spring_system.node_set;
  // boid_node_interact(boid_set,node_set,force_amp,force_k,vforce_amp);
  
  // set up some sliders
  repel_force_Tslider = new Tslider(0,'repel',0,repel_force); 
  flock.repel_force = repel_force;
  
  d_repel_Tslider = new Tslider(1, 'd_repel',0,d_repel);
  flock.d_repel = d_repel;
  
  align_force_Tslider = new Tslider(2, 'align',0,align_force);
  flock.align_force = align_force;
  
  d_align_Tslider = new Tslider(3, 'd_align',0,d_align);
  flock.d_align = d_align;
  
  gamma_node_Tslider = new Tslider(4, 'nodedamp',0,gamma_node);
  mass_spring_system.gamma_node = gamma_node;
  
  ks_Tslider = new Tslider(5, 'springk',ks/4,ks);
  mass_spring_system.ks = ks;
 
}

var dcount=0; // used for centroiding display

function draw() {

   let boid_set = flock.boid_set;
   let node_set = mass_spring_system.node_set;
   background(240);
    
   mass_spring_system.display_springs(); // display springs
   mass_spring_system.display_nodes();  // display nodes
   flock.display_flock();  // display boid flock
    
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
   
   // use slider info! update flock and springsystem params
   flock.repel_force = repel_force_Tslider.slider.value();
   flock.align_force = align_force_Tslider.slider.value();
   flock.d_repel = d_repel_Tslider.slider.value();
   flock.d_align = d_align_Tslider.slider.value();
   mass_spring_system.gamma_node = gamma_node_Tslider.slider.value();
   mass_spring_system.ks = ks_Tslider.slider.value();
   mass_spring_system.updateks(); // propagate to all springs
   repel_force_Tslider.text();
   align_force_Tslider.text();
   d_repel_Tslider.text();
   d_align_Tslider.text();
   gamma_node_Tslider.text();
   ks_Tslider.text();
   
}

// exponential short range forces between boids and nodes
// U propto force_amp*exp(-force_k*d) where d = distance between
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
      //let v = dv.mag;
      if (d < 3/force_k){
        let drhat = dr.copy();
        drhat.normalize(); //  is a unit vector
        let Force = drhat.copy(); // force in direction between
        Force.mult(-force_amp*exp(-force_k*d));
        if (vforce_amp >0){
          let dv = p5.Vector.sub(boidj.velocity,nodei.velocity);
          let vForce = dv.copy();
          vForce.mult(vforce_amp);// damping force depends on velocity difference
          Force.add(vForce);
        }
        let a_boid = Force.copy();  // equal and opposite forces
        let a_node = Force.copy();
        a_boid.div(-boidj.m); // acceleration on boid
        a_node.div( nodei.m); // acceleration on node
        boidj.acceleration.add(a_boid); // apply to accelerations
        nodei.acceleration.add(a_node);
      }
    }
  }
  
}
