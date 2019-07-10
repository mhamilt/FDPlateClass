//==============================================================================
// FD-Plate.js
//==============================================================================
const BoundaryCondition = {
  SIMPLE: 'simple',
  CLAMPED: 'clamped'
}
class FDPlate
{
  /** Constructor: Initialise with sample rate and FDPlate::PlateParameter struct
   or with default settings specifying nothing or just the sample
   rate. The Plate can be re-set using the setup() method*/
  constructor()
  {
    this.param = {
      /** Young's modulus*/
      youngs: 11e9,
      /** density (kg/m^3)*/
      density: 480,
      /** Poisson Ratios (< .5)*/
      poisson: .5,
      /** thickness (m)*/
      thickness: .003,
      /** x-axis plate length (m)*/
      lengthX: 1,
      /** y-axis plate length (m)*/
      lengthY: 1,
      /** T60 decay*/
      t60: 0.5,
      /** high frequency: percent of T60 (0 < tone < 1)*/
      tone: 0.2,
      /** boundary condtions*/
      bc: BoundaryCondition.SIMPLE,
      /**sample rate*/
      SR: 44100
    }
    this.state = {
      u: [],
      u1: [],
      u2: [],
      temp: []
    }

    //==========================================================================
    this.const = {
      pi: 3.14159265358979323846,
      maxGridSize: 3000,
      maxXgrid: 30.
    }

    this.coef = {
      /**x-axis plate length (m)*/
      Lx: null,
      /**y-axis plate length (m)*/
      Ly: null,
      /**readout position as percentage.*/
      rp: [1, 2, 3, 4],
      /**Excitation*/
      u0: null,
      v0: null,
      wid: null,
      ctr: [0.5, 0.5],

      //==========================================================================
      /**Loss coefficients*/
      sigma0: null,
      sigma1: null,
      //==========================================================================
      /**Scheme Coefficient*/
      A00: null,
      B00: null,
      B01: null,
      B11: null,
      B02: null,
      BC1: null,
      BC2: null,
      C00: null,
      C01: null,
      d0: null,
      //==========================================================================
      /**Derived Parameters*/
      Nx: null,
      Ny: null,
      ss: null,
      li: null,
      lo: null,
      lol: null,
      lor: null,
      /**Derived Parameters*/
      kappa: null,
      hmin: null,
      h: null,
      mu: null,
      k: null,
    }

  }

  //============================================================================
  //==========================================================================
  // Methods
  //==========================================================================
  signum(d)
  {
    return (d <= 0) ? 0 : 1;
  }
  /**
   set the plate to an initial condition of a raised cosine. This will overwrite
   all current values held on the plate.
   */
  addRc(ct)
  {
    if (arguments.length === 0)
    {
      ct = [0.5, 0.5];
    }

    {
      for (var xi = 1; xi < this.coef.Nx - 1; xi++)
      {
        let X = xi * h;

        for (var yi = 1; yi < this.coef.Ny - 1; yi++)
        {
          let cp = yi + (xi * this.coef.Ny);
          let Y = yi * h;
          let dist = Math.sqrt(Math.pow(X - (ct[0] * this.coef.Lx), 2) + Math.pow(Y - (ct[1] * this.coef.Ly), 2));
          let ind = signum((this.coef.wid * 0.5) - dist); // displacement (Math.logical)
          let rc = .5 * ind * (1 + Math.cos(2 * pi * dist / this.coef.wid)); // displacement
          this.state.u1[cp] += rc * this.coef.k;
        }
      }
    }
  }
  //==========================================================================
  /**
   Set the read-out point as a normalised position between 0 and 1.
   This position will be rounded to the nearest grid point

   @param xcoord x-axis co-ordinate
   @param ycoord y-axis co-ordinate
   */
  //  setOutputPosition( xcoord, ycoord);
  /**
   set the read-out position for interpolated output

   @param xCoord x-axis co-ordinate
   @param yCoord y-axis co-ordinate
   */
  //  setInterpOut(xCoord, yCoord);
  /**
   sets which function the is used when getting output
   */
  //  setOutputFunction(OutputMethod outType);
  //==========================================================================
  /**
   Get the output of the plate using the outputFunction pointer.

   @return returns value of which ever function outputFunction points to
   */
  // functio getOutput();

  /**
   Process an audio signal using the scheme as a plate reverb unit

   @param force the audio signal
   @return output from the plate at the desired point
   */
  reverb(force)
  {
    updateScheme();
    addForce(force);
    return getOutput();
  }
  //==========================================================================
  /**
   Update the time state of the scheme
   */
  updateScheme()
  {

    for (let xi = 2; xi < (this.coef.Nx - 2); ++xi)
    {
      for (let yi = 2; yi < (this.coef.Ny - 2); ++yi)
      {

        let cp = yi + (xi * this.coef.Ny);
        this.state.u[cp] = this.coef.B00 * this.state.u1[cp] +
          this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
          this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp + 2] + this.state.u1[cp - (2 * this.coef.Ny)] + this.state.u1[cp + (2 * this.coef.Ny)]) +
          this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
          this.coef.C00 * this.state.u2[cp] +
          this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
      }
    }

    // Update Side Boundaries
    //X-Axis
    for (let xi = 2; xi < this.coef.Nx - 2; ++xi)
    {
      //North
      {

        let cp = 1 + (xi * this.coef.Ny);
        this.state.u[cp] = this.coef.BC1 * this.state.u1[cp] +
          this.coef.B01 * (this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
          this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp + 2] + this.state.u1[cp - (2 * this.coef.Ny)] + this.state.u1[cp + (2 * this.coef.Ny)]) +
          this.coef.B11 * (this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny]) +
          this.coef.C00 * this.state.u2[cp] +
          this.coef.C01 * (this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
      }
      {
        //South

        let cp = this.coef.Ny - 2 + (xi * this.coef.Ny);
        this.state.u[cp] = this.coef.BC1 * this.state.u1[cp] +
          this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
          this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp - (2 * this.coef.Ny)] + this.state.u1[cp + (2 * this.coef.Ny)]) +
          this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
          this.coef.C00 * this.state.u2[cp] +
          this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
      }
    }

    // Y-Axis
    for (let yi = 2; yi < this.coef.Ny - 2; ++yi)
    {
      //West
      {
        let cp = yi + this.coef.Ny;
        this.state.u[cp] = this.coef.BC1 * this.state.u1[cp] +
          this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp + this.coef.Ny]) +
          this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp + 2] + this.state.u1[cp + (2 * this.coef.Ny)]) +
          this.coef.B11 * (this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
          this.coef.C00 * this.state.u2[cp] +
          this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp + this.coef.Ny]);
      }

      //East
      {
        let cp = yi + this.coef.Ny * (this.coef.Nx - 2);
        this.state.u[cp] = this.coef.BC1 * this.state.u1[cp] +
          this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny]) +
          this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp + 2] + this.state.u1[cp - (2 * this.coef.Ny)]) +
          this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny]) +
          this.coef.C00 * this.state.u2[cp] +
          this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny]);
      }
    }

    // Corner Boundaries
    {

      let cp = this.coef.Ny + 1;
      this.state.u[cp] = this.coef.BC2 * this.state.u1[cp] +
        this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
        this.coef.B02 * (this.state.u1[cp + 2] + this.state.u1[cp + (2 * this.coef.Ny)]) +
        this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
        this.coef.C00 * this.state.u2[cp] +
        this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
    }
    {
      let cp = 2 * (this.coef.Ny - 1);
      this.state.u[cp] = this.coef.BC2 * this.state.u1[cp] +
        this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
        this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp + (2 * this.coef.Ny)]) +
        this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
        this.coef.C00 * this.state.u2[cp] +
        this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
    }
    {
      let cp = this.coef.Ny * (this.coef.Nx - 2) + 1;
      this.state.u[cp] = this.coef.BC2 * this.state.u1[cp] +
        this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
        this.coef.B02 * (this.state.u1[cp + 2] + this.state.u1[cp - (2 * this.coef.Ny)]) +
        this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
        this.coef.C00 * this.state.u2[cp] +
        this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
    }
    {
      let cp = this.coef.Ny * (this.coef.Nx - 1) - 2;
      this.state.u[cp] = this.coef.BC2 * this.state.u1[cp] +
        this.coef.B01 * (this.state.u1[cp - 1] + this.state.u1[cp + 1] + this.state.u1[cp - this.coef.Ny] + this.state.u1[cp + this.coef.Ny]) +
        this.coef.B02 * (this.state.u1[cp - 2] + this.state.u1[cp - (2 * this.coef.Ny)]) +
        this.coef.B11 * (this.state.u1[cp - 1 - this.coef.Ny] + this.state.u1[cp + 1 - this.coef.Ny] + this.state.u1[cp + 1 + this.coef.Ny] + this.state.u1[cp - 1 + this.coef.Ny]) +
        this.coef.C00 * this.state.u2[cp] +
        this.coef.C01 * (this.state.u2[cp - 1] + this.state.u2[cp + 1] + this.state.u2[cp - this.coef.Ny] + this.state.u2[cp + this.coef.Ny]);
    }

    // swap pointers
    let dummyptr = this.state.u2;
    this.state.u2 = this.state.u1;
    this.state.u1 = this.state.u;
    this.state.u = dummyptr;

  }
  /**
   add force to the relevant section of the plate

   @param force force in Newtons
   */
  addForce(force)
  {
    this.state.u1[this.coef.li] += this.coef.d0 * force;
  }
  //==========================================================================
  /**
   Print plate parameter information
   */
  //  printInfo();
  /**
   Print Internal Coefficients
   */
  //  printCoefs();

  //==========================================================================
  /**
   setup the plate with a given sample rate and boundary condition type: still under construction

   @param sampRate Sample Rate in Hz
   */
  setup(sampRate)
  {
    // I/O Parameters
    this.coef.rp[0] = .5;
    this.coef.rp[1] = .4;
    this.coef.rp[2] = .5;
    this.coef.rp[3] = .5; // readout position as percentage.

    //Excitation
    this.coef.ctr[0] = .5;
    this.coef.ctr[1] = .5; // centre point of excitation as percentage
    this.coef.wid = .25; // width (m)
    this.state.u0 = 0;
    this.coef.v0 = 1; // excitation displacement and velocity

    let E = this.param.youngs;
    let H = this.param.thickness;
    let nu = this.param.poisson;
    let rho = this.param.density;
    let D = (E * Math.pow(H, 3)) / (12 * (1 - Math.pow(nu, 2)));
    this.coef.kappa = Math.sqrt(D / (rho * H));
    this.coef.SR = sampRate; // internal class sampling rate
    this.coef.k = 1 / this.coef.SR; // time step

    this.coef.Lx = this.param.lengthX;
    this.coef.Ly = this.param.lengthY;
    this.setLoss(this.param.t60, this.param.tone);
    this.setGridSpacing();
    this.setCoefs(this.param.bcType, rho, H);

    this.coef.d0 = (Math.pow(this.coef.k, 2)) / (rho * H * (Math.pow(this.coef.h, 2))) * (1 / (1 + this.coef.k * this.coef.sigma0))

    //Set Input and Output Indeces
    this.coef.li = Math.floor((this.coef.Ny * (this.coef.ctr[1] * this.coef.Nx)) + (this.coef.ctr[0] * this.coef.Ny));
    this.coef.lo = Math.floor((this.coef.Ny * (this.coef.rp[1] * this.coef.Nx)) + (this.coef.rp[0] * this.coef.Ny));

    //    Update flags

    this.state.u = new Array(this.coef.ss);
    this.state.u1 = new Array(this.coef.ss);
    this.state.u2 = new Array(this.coef.ss);
    for (let j = 0; j < this.state.u.length; j++)
    {
      this.state.u[j] = 0;
      this.state.u1[j] = 0;
      this.state.u2[j] = 0;
    }
  }
  //==========================================================================
  /**
   get the velocity output from the plate. Rounding to the nearest grid point

   @return velocity at specified read-out point
   */
  //  getVelocityOutput();
  /**
   Gets the amplitude output from the plate. Rounding to the nearest grid point

   @return amplitude at specified read-out point
   */
  //  getAmplitudeOutput();
  /**
   Use the internal interpolation method to calculate the output from the plate

   @return interpolated value at specified read-out point
   */
  //  getInterpOut();
  //==========================================================================
  /**
   Populates the internal interpolation lookup table.
   */
  //  setInterpLookTable();
  //==========================================================================
  /**
   Sets the loss coefficients

   @param t60 the T60 decay in seconds
   @param tone percentage of high frequency decay relative to T60, between .1 and 1.
   Values outside this range will be capped
   */
  setLoss(t60, tone)
  {
    tone = Math.min(Math.max(tone, 0.1), 0.99);
    let lowFrequencyBand = 100;
    let highFrequencyBand = 1000;
    let highT60 = t60 * tone;
    let z1 = 2 * this.coef.kappa * (2 * this.const.pi * lowFrequencyBand) / (2 * Math.pow(this.coef.kappa, 2));
    let z2 = 2 * this.coef.kappa * (2 * this.const.pi * highFrequencyBand) / (2 * Math.pow(this.coef.kappa, 2));
    this.coef.sigma0 = 6 * Math.log(10) * (-z2 / t60 + z1 / highT60) / (z1 - z2);
    this.coef.sigma1 = 6 * Math.log(10) * (1 / t60 - 1 / highT60) / (z1 - z2);
  }
  /**
   Sets the coeffiecients for the scheme

   @param bcType Boundary Condition
   @param rho plate material density
   @param H plate thickness
   */
  setCoefs(bcType, rho, H)
  {

    // currentBoundCon = bcType; //update flag

    // coefficients are named based on position on the x and y axes.
    this.coef.A00 = 1 / (1 + this.coef.k * this.coef.sigma0); // Central Loss Coeffient (INVERTED)

    //// Current time step (B) coeffients
    // There are six unique coefficients for B coefs
    this.coef.B00 = (-Math.pow(this.coef.mu, 2) * 20 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 + 2) * this.coef.A00; // center
    this.coef.B01 = (-Math.pow(this.coef.mu, 2) * -8 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2))) * this.coef.A00; // 1-off
    this.coef.B11 = (-Math.pow(this.coef.mu, 2) * 2) * this.coef.A00; // diag
    this.coef.B02 = (-Math.pow(this.coef.mu, 2) * 1) * this.coef.A00; // 2-off

    switch (bcType)
    {
      case BoundaryCondition.CLAMPED:
      {
        this.coef.BC1 = (-Math.pow(this.coef.mu, 2) * 21 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 + 2) * this.coef.A00; // Side
        this.coef.BC2 = (-Math.pow(this.coef.mu, 2) * 22 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 + 2) * this.coef.A00; // Corner
        break;
      }
      case BoundaryCondition.SIMPLE:
      default:
      {
        this.coef.BC1 = (-Math.pow(this.coef.mu, 2) * 19 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 + 2) * this.coef.A00; // Side
        this.coef.BC2 = (-Math.pow(this.coef.mu, 2) * 18 + (2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 + 2) * this.coef.A00; // Corner
        break;
      }
    }

    // Previous time step (C) coeffients
    this.coef.C00 = (-(2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * -4 - (1 - this.coef.sigma0 * this.coef.k)) * this.coef.A00;
    this.coef.C01 = -(2 * this.coef.sigma1 * this.coef.k / Math.pow(this.coef.h, 2)) * this.coef.A00;

    //input force coefficient
    this.coef.d0 = Math.pow(this.coef.k, 2) / (this.coef.rho * this.coef.H * Math.pow(this.coef.h, 2)) * (1 / (1 + this.coef.k * this.coef.sigma0)) * this.coef.A00;

  }
  /**
   Sets the grid spacing for the scheme
   */
  setGridSpacing()
  {
    // stability condition
    this.coef.hmin = (Math.sqrt(4 * this.coef.k * (this.coef.sigma1 + Math.sqrt(Math.pow(this.coef.sigma1, 2) + Math.pow(this.coef.kappa, 2)))));
    this.coef.Nx = Math.floor(this.coef.Lx / this.coef.hmin); // number of segments x-axis
    this.coef.Ny = Math.floor(this.coef.Ly / this.coef.hmin); // number of segments y-axis

    this.coef.h = Math.sqrt(this.coef.Lx * this.coef.Ly / (this.coef.Nx * this.coef.Ny));; // adjusted grid spacing x/y

    this.coef.Nx = this.coef.Nx + 1;
    this.coef.Ny = this.coef.Ny + 1; // grid point number x and y
    this.coef.mu = (this.coef.kappa * this.coef.k) / Math.pow(this.coef.h, 2); // scheme parameter

    this.coef.ss = this.coef.Nx * this.coef.Ny; // total grid size.
  }
  //==========================================================================
  printvars()
  {
    console.log("--- Coefficient Info --- \n\n");
    console.log("Loss A		: %.4f \n", this.coef.A00);
    console.log("Centre B    : %.4f \n", this.coef.B00);
    console.log("1-Grid B    : %.4f \n", this.coef.B01);
    console.log("2-Grid B	: %.4f \n", this.coef.B02);
    console.log("Diagonal B  : %.4f \n", this.coef.B11);
    console.log("Centre C	: %.4f \n", this.coef.C00);
    console.log("1-Grid C    : %.4f \n", this.coef.C01);
    console.log("Side Bound	: %.4f \n", this.coef.BC1);
    console.log("Cornr Bound : %.4f \n", this.coef.BC2);

    console.log("\n--- Scheme Info --- \n\n");
    console.log("Size		: %.1fm2 \n", this.coef.Nx * this.coef.h * this.coef.Ny * this.coef.h);
    console.log("Grid X-Ax   : %d \n", this.coef.Nx);
    console.log("Grid Y-Ax   : %d \n", this.coef.Ny);
    console.log("Total P		: %d \n", this.coef.ss);
    console.log("Dur(samps)	: %d \n", this.coef.Nf);
    console.log("In_cell		: %d\n", this.coef.li);
    console.log("Out_cell	: %d\n", this.coef.lo);
    console.log("Youngs		: %.2e\n", this.coef.E);
    console.log("Sigma 0		: %f\n", this.coef.sigma0);
    console.log("Sigma 1		: %f\n", this.coef.sigma1);
  }
}
