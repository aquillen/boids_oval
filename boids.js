function Boid(dxx) {
  // constructor
  this.acceleration = createVector(0,0);
  let vw = 1;
  this.velocity = createVector(random(-vw,vw),random(-vw,vw));
  let xw = 130;
  this.position = createVector(random(-xw,xw),random(-xw,xw));
  this.r = 5.0;  // for display and in pixels
  this.m = boid_mass;
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
