

// create a slider and store a label and position for it so we
// can update text labeling during draw
// each slider should have a different index so they are 
// all in a line
function Tslider(index,label,min,value){
  const slider_len = 80;
  const slider_len_s = '80px';
  const x0 = 10;  
  const sep = 10; // separation between sliders
  const y0 = 10
  this.label = label;
  this.slider_x = index*(slider_len + sep) + x0;
  this.slider_y = y0;
  // this.slider_w = slider_len;
  let max = value*2;
  let step = (max-min)/20;
  
  // constructor
  this.slider = createSlider(min,max,value,step);
  this.slider.position(this.slider_x, this.slider_y);
  this.slider.style('width',  slider_len_s);
  
}

// label the slider
Tslider.prototype.text = function(){
   fill(0); // black
   text(this.label,this.slider_x + 10,this.slider_y+20); // offset the label
}
