
function Flock(nboids) {
  // An array for all the boids
  this.boid_set = []; // Initialize the array of boids
  this.dx = dx;
  // this.maxforce = 20;
  // this.maxspeed = 2;
  this.boid_mass = boid_mass; // boid mass for flock
    
  this.boidspeed = boidspeed;
  
  this.d_attract = d_attract; // distance for attraction/repulsion forces between Boids
  this.d_repel   = d_repel;  // distance for attraction/repulsion forces between Boid
  this.d_align   = d_align;  // distance for alignment
  
  this.attract_force = attract_force;  // attact force size
  this.repel_force   = repel_force;    // repel force size
  this.align_force   = align_force;
  this.propel_force  = propel_force;

  
  // create nboids of Boids!
  for (let i=0;i<nboids;i++){
    let b = new Boid(this.dx,this.boid_mass);
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
   if (this.repel_force ==0){
      return;
   }
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
     if (this.attract_force ==0){
       return;
     }
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
           if (count >0){
              pos_sum.div(count);   // is average of nearby boid positions!
              let Force = p5.Vector.sub(pos_sum, bi.position); // target direction
              Force.normalize();
              Force.mult(this.boidspeed); // boidspeed used here!
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
     if (this.align_force==0){
       return;
     }
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
          v_ave.mult(this.boidspeed); // boidspeed used here!
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
     // boid_node_interact(boid_set,node_set,force_amp,force_k,vforce_amp); 
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



function Boid(dxx,boid_mass) {
  // constructor
  this.acceleration = createVector(0,0);
  const vw = 1;
  this.velocity = createVector(random(-vw,vw),random(-vw,vw));
  const xw = 130;
  this.position = createVector(random(-xw,xw),random(-xw,xw));
  this.r = 5.0;  // for display and in pixels
  this.m = boid_mass;  // global!
  this.dx = dxx; // scale multiply by this to get pixels

}

// Lots of Boid stuff from p5 examples
Boid.prototype.render = function() {
  // Draw a triangle rotated in the direction of velocity
  var theta = this.velocity.heading() + radians(90);
  push();
  fill(0,0,100);
  stroke(20);
  strokeWeight(1);
  translate(0,0);
  translate(width/2+this.position.x*this.dx,height/2+this.position.y*this.dx);
  rotate(theta);
  beginShape();
  vertex(0, -this.r*2);
  vertex(-this.r, this.r*2);
  vertex(this.r, this.r*2);
  endShape(CLOSE);
  pop();

}


// wrap boids, but note that forces are not yet wrapped.
Boid.prototype.borders = function(){
  //print(this.dx);
  let wreal = (width + this.r)/this.dx ;
  let hreal = (height+ this.r)/this.dx ;
  if (this.position.x > 0.5*wreal){
    this.position.x -= wreal;
  }
  if (this.position.x < -0.5*hreal){
    this.position.x += wreal;
  }
  if (this.position.y > 0.5*wreal){
    this.position.y -= hreal;
  }
  if (this.position.y < -0.5*hreal){
    this.position.y += hreal;
  }
}
