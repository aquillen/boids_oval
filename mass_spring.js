


// set up a mass spring system
// depends on ks, gammas, nnodes, big_radius, dt, mass
function Mass_spring_system(nnodes,big_radius,node_mass){
  
  // constructor
  this.node_set = []; // Initialize the node array
  this.spring_set = []; // Initialize the spring array
  this.nsprings = 0;  // number of springs
  this.eps = 1e-6; // smoothing length in constructor ..
  this.ks = ks; // to set spring constants
  this.gammas = gammas; // to set dissipation rates
  this.gamma_node = gamma_node;
  this.dx = dx;
  this.node_mass = node_mass;

  // construct nodes in a circle, node origin is center of canvas
  for (let i = 0; i < nnodes; i++) {
    let b = new Node(big_radius*cos(i*2.0*PI/nnodes),
                     big_radius*sin(i*2.0*PI/nnodes),this.dx);
    b.m = this.node_mass; // node mass
    this.node_set.push(b);
  }
  
  // construct springs with nearest neighbor separation integer nene
  // ks are adjusted with ksfac*ks
  // damping is adjusted with gammafac*gammas
  this.add_spring_set = function(nene,ks0,gammas0){
    for (let k = 0; k < nnodes; k++) {  
      let s = new Spring(k,(k+nene)%nnodes,ks0,gammas0);  // which nodes the spring connects
      this.spring_set.push(s);
      this.nsprings++;
    }
  }
  // in constructor, three sets of springs added
  this.add_spring_set(1,1.0*this.ks,1*this.gammas); // add nearest neighbor springs
  this.add_spring_set(2,1.0*this.ks,1*this.gammas); // add springs to every other node
  this.add_spring_set(3,1.0*this.ks,1*this.gammas); // add springs to every 3rd node
  
  // set spring k to its rest length
  this.rest_spring = function(k){
    let i = this.spring_set[k].i;
    let j = this.spring_set[k].j;
    let nodei = this.node_set[i];
    let nodej = this.node_set[j];
    let distij = p5.Vector.dist(nodei.position,nodej.position);
    this.spring_set[k].L0 = distij;
  }
    
  // set rest spring lengths to current lengths for all springs
  this.set_spring_lengths = function(){
    for(let k=0;k< this.nsprings; k++){
      this.rest_spring(k);
    }
  }
  
  // if you want to change ks for all springs
  this.updateks = function(){
    for(let k=0;k< this.nsprings; k++){
      this.spring_set[k].ks = this.ks;  // note assumes all springs the same ks!
    }
  }
      
  // for adding a small wedge to the spring system
  this.add_wedge = function(i0,di){
      let n = this.node_set.length;
      let indexi = i0;
      let nodei = this.node_set[indexi];
      let indexj = indexi+di;
      let nodej = this.node_set[indexj];
      let mid = p5.Vector.add(nodei.position,nodej.position);
      mid.normalize();  // point direction in between 
      let ri = nodei.position.mag();
      let dist = p5.Vector.dist(nodei.position,nodej.position);
      mid.mult(ri-dist);
      let b = new Node(mid.x,mid.y,this.dx);
      b.m = nodei.m; // node mass
      this.node_set.push(b);
      let nodek = this.node_set[n]; // the new node
 
      let s = new Spring(indexi,n,this.ks,this.gammas); 
      this.spring_set.push(s);
      this.rest_spring(this.nsprings); 
      this.nsprings++;
      s = new Spring(indexj,n,this.ks,this.gammas); 
      this.spring_set.push(s);
      this.rest_spring(this.nsprings); 
      this.nsprings++;
      // this.set_spring_lengths();
    }
    
// display nodes
  this.display_nodes = function(){
    stroke(51);
    strokeWeight(1);
    //fill('red');
    fill('rgb(255,0,0)');
    let n = this.node_set.length; // number of nodes
    for (let i=0;i<n;i++){
      ellipse(width/2+this.node_set[i].position.x*this.dx, 
              height/2+this.node_set[i].position.y*this.dx,10,10);
    }
  }

// display springs
  this.display_springs = function(){
    stroke(100,100,0);
    strokeWeight(3);
    for (let i =0; i< this.nsprings; i++){
      let node_i = this.spring_set[i].i;
      let node_j = this.spring_set[i].j;
      // print(node_i,node_j);
      line(width/2 +this.node_set[node_i].position.x*this.dx,
        height/2+this.node_set[node_i].position.y*this.dx,
        width/2 +this.node_set[node_j].position.x*this.dx,
        height/2+this.node_set[node_j].position.y*this.dx); 
    }
  }

  // compute spring force for 1 spring on the 2 nodes it affects
  // k is the spring index
  // add to accelerations on the two nodes
  this.springforce = function(k){
    let i = this.spring_set[k].i;
    let j = this.spring_set[k].j;
    let nodei = this.node_set[i];
    let nodej = this.node_set[j];
    let ks = this.spring_set[k].ks;
    let gammas = this.spring_set[k].gammas;
    let L0 = this.spring_set[k].L0;
    
    let Lvec = p5.Vector.sub(nodei.position, nodej.position); // difference dx vector
    let Lhat = Lvec.copy();
    Lhat.normalize(); // unit vector!
    let L = Lvec.mag(); // length

    // kForce is spring force
    let kForce = p5.Vector.mult(Lhat,ks*(L-L0));  // is a vector in the direction between nodes
    let kai = p5.Vector.mult(kForce,1.0/nodei.m);  // acceleration on node i
    let kaj = p5.Vector.mult(kForce,1.0/nodej.m);  // acceleration on node j
    // sum spring accelerations on nodes
    nodei.acceleration.sub(kai); // store by adding to node acceleration
    nodej.acceleration.add(kaj);
    
    // dForce is damping force
    if (gammas>0){
      let dV = p5.Vector.sub(nodei.velocity, nodej.velocity);  // dV vector
      let dxdotdv = p5.Vector.dot(Lvec,dV); // dx dot dv
      let mu_mass = nodei.m * nodej.m/(nodei.m + nodej.m); // reduced mass
      let dForce = p5.Vector.mult(Lhat,gammas*dxdotdv*mu_mass/(L+ this.eps));
      let dai = p5.Vector.mult(dForce,1.0/nodei.m);  // acceleration on node i
      let daj = p5.Vector.mult(dForce,1.0/nodej.m);  // acceleration on node j
      // sum damping accelerations on nodes
      nodei.acceleration.sub(dai); // store by adding to node acceleration
      nodej.acceleration.add(daj);
    }
    
  }

  // compute spring forces on nodes for all springs, including damping
  this.applysprings = function(){
    for(let k = 0;k< this.nsprings;k++){
      this.springforce(k);  // adds to accelerations
    }
  }

  // zero all the accelerations on the nodes
  this.zeroaccel = function(){
    let n = this.node_set.length;
    for(let k = 0;k< n;k++){
      this.node_set[k].acceleration.set(0,0);
    }
  }
  
  this.centroid = function(){
    let sumvec = createVector(0,0);
    let n = this.node_set.length;
    for(let i = 0;i< n;i++){
      nodei = this.node_set[i];
      sumvec.add(nodei.position);
    }
    if (n>0){
      sumvec.div(n);
    }
    return sumvec;
  }
    
  this.shift = function(centroid){
    let n = this.node_set.length;
    for(let i = 0;i< n;i++){
            nodei = this.node_set[i];
            nodei.position.sub(centroid);
     }
  }
  
  // damping directly on nodes
  this.node_damp = function(){
    let n = this.node_set.length;
    for(let i = 0;i< n;i++){
      let nodei = this.node_set[i];
      let dForce = p5.Vector.mult(nodei.velocity,this.gamma_node);
      nodei.acceleration.sub(dForce);
    }
  }
  
  // compute all accelerations on nodes!
  this.compute_accel = function(boid_set,node_set,
                                 force_amp,force_k,vforce_amp){
    // this.zeroaccel(); // zero accelerations
    this.applysprings(); // compute accelerations from springs
    this.node_damp(); // apply damping direcly to nodes
    // boid_node_interact(boid_set,node_set,force_amp, force_k,vforce_amp);
  }

  this.single_timestep = function(dt){ 
    // accelerations have to be already computed
    let n = this.node_set.length;
    for (let i=0;i<n;i++){
      let nodei = this.node_set[i];
      let dv = p5.Vector.mult(nodei.acceleration,dt); 
      nodei.velocity.add(dv); // 
      let dr = p5.Vector.mult(nodei.velocity,dt);
      nodei.position.add(dr); // 
    }  
  }
  
  // part of constructor for mass spring system!
  this.set_spring_lengths();  // set rest lengths!
  // this.node_set[0].position.add(10,0); // perturb one node!
  this.compute_accel(); // compute accelerations from springs

}


// a single mass node
function Node(x,y,dxx){
  this.position = createVector(x,y);
  this.velocity = createVector(0,0);
  this.acceleration = createVector(0,0);
  this.m = node_mass;
  this.dx = dxx; // scale
}

// a single spring  // class here is equivalent to function above
function Spring(i,j,ks0,gammas0){
  this.i = i;  // spring connects mass node i to mass node j
  this.j = j;
  this.ks = ks0; // spring constant
  this.gammas = gammas0; // damping coeff
  this.L0 = 1.0;  // rest length
}
