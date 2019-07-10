// Main Sketch
let screensize = (4 * $("#jumbo-canvas").width()) / 5;
let scale;
let shdr;
let zscl = 1.5e5;
//------------------------------------------------------------------------------
let plate = new FDPlate();
plate.setup(44.1e3);
//------------------------------------------------------------------------------
// Display
let display_3d = false;
let frameSkip = 10;
//------------------------------------------------------------------------------
// Input
let input_sine = true;
let radpersec;
let cur_rad = 0
let forceAmp = 10;
let forceIn = {
  x: 0,
  y: 0
};
let input_rect_timer = 0
//------------------------------------------------------------------------------
preload = function()
{
  plateShader = loadShader('js/shaders/vert.shader', 'js/shaders/frag.shader');
};
//------------------------------------------------------------------------------
function setup()
{
  //----------------------------------------------------------------------------
  setForce(100)
  //----------------------------------------------------------------------------
  var canvas = createCanvas(screensize, screensize, WEBGL);
  frameRate(30);
  // gl = canvas.getContext('webgl');
  //----------------------------------------------------------------------------
  canvas.parent('sketch-holder');
  pg = createGraphics(200, 200);
  pg.textSize(75);
  pg.background(0, 100);
  fill(255);
  texture(pg);
  shader(plateShader);
  plateShader.setUniform('zscale', 0.0)
  //----------------------------------------------------------------------------
  forceIn = {
    x: plate.coef.ctr[0] * width,
    y: plate.coef.ctr[1] * height
  }
}
//------------------------------------------------------------------------------
function draw()
{
  drawPlate();
}
//------------------------------------------------------------------------------
function windowResized()
{
  screensize = (4 * $("#jumbo-canvas").width()) / 5;
  resizeCanvas(screensize, screensize);
}
//------------------------------------------------------------------------------
function drawPlate()
{
  //----------------------------------------------------------------------------
  for (var i = 0; i < frameSkip; ++i)
  {
    if (input_sine)
    {
      plate.addForce(sin(cur_rad) * forceAmp)
      cur_rad += radpersec;
    }
    plate.updateScheme();
  }
  //----------------------------------------------------------------------------
  background(50);
  scale = width / (plate.coef.Nx - 1);
  push();
  shader(plateShader);
  //----------------------------------------------------------------------------
  if (display_3d)
  {
    plateShader.setUniform('zscale', 1.0)
    translate(0, 0, -200);
    rotateX(PI / 3);
    rotateZ(millis() * 5e-4);
    rotateY(millis() * 5e-4);
  }
  else
  {
    plateShader.setUniform('zscale', 0.0)
  }
  translate(-width / 2, -height / 2, 0);
  for (var yi = 0; yi < plate.coef.Ny; ++yi)
  {
    beginShape(TRIANGLE_STRIP);
    for (var xi = 0; xi < plate.coef.Nx; ++xi)
    {
      var cp = (xi) + ((yi) * plate.coef.Nx); // current povar
      vertex(xi * scale, yi * scale, Math.floor(plate.state.u[cp] * zscl));
      vertex(xi * scale, (yi + 1) * scale, Math.floor(plate.state.u[cp + plate.coef.Ny] * zscl));
    }
    endShape()
  }
  //----------------------------------------------------------------------------
  rectMode(CENTER)
  let col = 70 + (sin(cur_rad) + 1) * 50;
  if (!display_3d && input_sine)
  {
    if (millis() - input_rect_timer < 5000)
    {
      fill(0, col, col)
      rect(forceIn.x, forceIn.y, 20, 20)
    }
  }
  //----------------------------------------------------------------------------
  pop();
  //----------------------------------------------------------------------------

}
//------------------------------------------------------------------------------
function mousePressed()
{
  if (
    mouseX > 0 &&
    mouseX < width &&
    mouseY > 0 &&
    mouseY < height
  )
  {
    input_rect_timer = millis();
    var point = constrain(Math.floor(((plate.coef.Ny - 1) * mouseX / height)), 1, plate.coef.Ny - 2) + plate.coef.Ny
    var c = [
      mouseY / height,
      mouseX / width
    ]
    forceIn = {
      x: mouseX,
      y: mouseY
    }
    if (!input_sine)
    {
      plate.addRc(c)
    }
  }
}
//------------------------------------------------------------------------------
function mouseDragged()
{
  if (
    mouseX > 0 &&
    mouseX < width &&
    mouseY > 0 &&
    mouseY < height
  )
  {
    input_rect_timer = millis();
    var point = constrain(Math.floor(((plate.coef.Ny - 1) * mouseX / height)), 1, plate.coef.Ny - 2) + plate.coef.Ny
    forceIn = {
      x: mouseX,
      y: mouseY
    }
    plate.coef.li = point;
  }
}
//------------------------------------------------------------------------------
function setForce(freq)
{
  frameSkip = (freq > 500) ? 1 : 10
  radpersec = TAU * freq / plate.coef.SR
}
