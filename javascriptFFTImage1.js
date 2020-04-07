<script>
var Fourier = (function () {
      /******************
       * work functions */
      function filter(data, dims, lowPass, highPass) {
        var lowPassSq = Math.pow(lowPass, 2);
        var highPassSq = Math.pow(highPass, 2);
        var N = dims[1];
        var M = dims[0];
        for (var k = 0; k < N; k++) {
          for (var l = 0; l < M; l++) {
            var idx = k*M + l;
            var d = Math.pow(k-M/2, 2) + Math.pow(l-N/2, 2);
            if (
              d > lowPassSq && isNaN(highPass) ||
              d < highPassSq && isNaN(lowPass) ||
              d < lowPassSq && !isNaN(lowPass) && !isNaN(highPass) ||
              d > highPassSq && !isNaN(lowPass) && !isNaN(highPass)
            ) {
              data[idx] = new Fourier.Complex(0, 0);
            }
          }
        }
      }

      function FFT(sig, out) {
        rec_FFT(out, 0, sig, 0, sig.length, 1);
      }

      function rec_FFT(out, start, sig, offset, N, s) {
        if (N === 1) {
          out[start] = new Complex(sig[offset], 0); // array
        } else {
          rec_FFT(out, start, sig, offset, N/2, 2*s);
          rec_FFT(out, start+N/2, sig, offset+s, N/2, 2*s);
          for (var k = 0; k < N/2; k++) {
            var twiddle = cisExp(-2*Math.PI*k/N);
            var t = out[start+k];
            out[start+k] = t.plus(twiddle.times(out[start+k+N/2]));
            out[start+k+N/2] = t.minus(
              twiddle.times(out[start+k+N/2])
            );
          }
        }
      }

      function invFFT(transform, sig) {
        rec_invFFT(sig, 0, transform, 0, transform.length, 1);
        for (var ai = 0; ai < sig.length; ai++) {
          sig[ai] = sig[ai].real/sig.length;
        }
      }

      function rec_invFFT(sig, start, transform, offset, N, s) {
        if (N === 1) {
          sig[start] = transform[offset];
        } else {
          rec_invFFT(sig, start, transform, offset, N/2, 2*s);
          rec_invFFT(sig, start+N/2, transform, offset+s, N/2, 2*s);
          for (var k = 0; k < N/2; k++) {
            var twiddle = cisExp(2*Math.PI*k/N);
            var t = sig[start+k];
            sig[start+k] = t.plus(twiddle.times(sig[start+k+N/2]));
            sig[start+k+N/2] = t.minus(
              twiddle.times(sig[start+k+N/2])
            );
          }
        }
      }

      /********************
       * helper functions */
      function cisExp(x) { // e^ix = cos x + i*sin x
        return new Complex(Math.cos(x), Math.sin(x));
      }

      /***********
       * objects */
      function Complex(re, im) {
        this.real = re;
        this.imag = im;
      }
      Complex.prototype.magnitude2 = function() {
        return this.real*this.real + this.imag*this.imag;
      };
      Complex.prototype.magnitude = function() {
        return Math.sqrt(this.magnitude2());
      };
      Complex.prototype.plus = function(z) {
        return new Complex(this.real+z.real, this.imag+z.imag);
      };
      Complex.prototype.minus = function(z) {
        return new Complex(this.real-z.real, this.imag-z.imag);
      };
      Complex.prototype.times = function(z) {
        if (typeof z === 'object') { // complex multiplication
          var rePart = this.real*z.real - this.imag*z.imag;
          var imPart = this.real*z.imag + this.imag*z.real;
          return new Complex(rePart, imPart);
        } else { // scalar multiplication
          return new Complex(z*this.real, z*this.imag);
        }
      };

      return {
        Complex: Complex,
        transform: FFT,
        invert: invFFT,
        filter: filter
      };
    })();

var FourierImageAnalysis = (function() {
  /**********
   * config */
  var dims = [-1, -1]; // will be set later
  var cc = 9e-3; // contrast constant

  /*********************
   * working variables */
  var canvases;
  var ctxs;
  var h;
  var $h; // h hat
  var h_; // h prime, the reconstructed h values

  /******************
   * work functions */
  function initFourierImageAnalysis() {
    // event listeners
    $s('#draw-cs-btn').addEventListener('click', function() {
      loadImage('cs.png');
    });

    $s('#draw-circle-btn').addEventListener('click', function() {
      loadImage('circle.png');
    });

    $s('#draw-grace-btn').addEventListener('click', function() {
      loadImage('grace.png');
    });

    $s('#draw-img-btn').addEventListener('click', function() {
      loadImage($s('#img-url').value);
    });

    $s('#transform-btn').addEventListener('click', function() {
      var start = +new Date();

      if (!h()) {
        return alert(
          'You need to draw an image to canvas 1 first.'
        );
      }

      // placed in a callback so the UI has a chance to update
      disableButtons(transformAction);
    });

    $s('#reconstruct-btn').addEventListener('click', function() {
      var start = +new Date();

      if (!$h()) {
        return alert(
          'You first need to compute the Fourier transform.'
        );
      }

      // placed in a callback so the UI has a chance to update
      disableButtons(reconstructAction);
    });

    $s('#difference-btn').addEventListener('click', function() {
      var start = +new Date();

      if (!h_()) {
        return alert('You haven\'t reconstructed an image yet.');
      }

      // placed in a callback so the UI has a chance to update
      disableButtons(differenceAction);
    });

    // initialize the working variables
    canvases = [], ctxs = [];
    h = $h = h_ = function() { return false; };
  }

  function loadImage(loc) {
    var start = +new Date();

    // placed in a callback so the UI has a chance to update
    disableButtons(function() {
      // draw the initial image
      var img = new Image();
      img.addEventListener('load', function() {
        // make each canvas the image's exact size
        dims[0] = img.width;
        dims[1] = img.height;
        for (var ai = 0; ai < 4; ai++) {
          canvases[ai] = $s('#canvas'+ai);
          canvases[ai].width = dims[0];
          canvases[ai].height = dims[1];
          ctxs[ai] = canvases[ai].getContext('2d');
        }

        // draw the image to the canvas
        ctxs[0].drawImage(img, 0, 0, img.width, img.height);

        // grab the pixels
        var imageData = ctxs[0].getImageData(
          0, 0, dims[0], dims[1]
        );
        var h_es = []; // the h values
        for (var ai = 0; ai < imageData.data.length; ai+=4) {
          // greyscale, so you only need every 4th value
          h_es.push(imageData.data[ai]);
        }

        // initialize the h values
        h = function(n, m) {
          if (arguments.length === 0) return h_es;

          var idx = n*dims[0] + m;
          return h_es[idx];
        }; // make it a function so the code matches the math

        enableButtons();

        var duration = +new Date() - start;
        console.log(
          'It took ' + duration + 'ms to draw the image.'
        );
      });
      img.crossOrigin = "anonymous";
      img.src = loc;
    });
  }

  function transformAction() {
    // compute the h hat values
    var h_hats = [];
    Fourier.transform(h(), h_hats);
    h_hats = Fourier.shift(h_hats, dims);

    // get the largest magnitude
    var maxMagnitude = 0;
    for (var ai = 0; ai < h_hats.length; ai++) {
      var mag = h_hats[ai].magnitude();
      if (mag > maxMagnitude) {
        maxMagnitude = mag;
      }
    }

    // apply a low or high pass filter
    var lowPassRadius = parseInt(
      $s('#low-freq-radius').value
    ); // low pass radius
    var highPassRadius= parseInt(
      $s('#high-freq-radius').value
    ); // high pass radius
    Fourier.filter(h_hats, dims, lowPassRadius, highPassRadius);

    // store them in a nice function to match the math
    $h = function(k, l) {
      if (arguments.length === 0) return h_hats;

      var idx = k*dims[0] + l;
      return h_hats[idx];
    };

    // draw the pixels
    var currImageData = ctxs[1].getImageData(
      0, 0, dims[0], dims[1]
    );
    var logOfMaxMag = Math.log(cc*maxMagnitude+1);
    for (var k = 0; k < dims[1]; k++) {
      for (var l = 0; l < dims[0]; l++) {
        var idxInPixels = 4*(dims[0]*k + l);
        currImageData.data[idxInPixels+3] = 255; // full alpha
        var color = Math.log(cc*$h(l, k).magnitude()+1);
        color = Math.round(255*(color/logOfMaxMag));
        // RGB are the same -> gray
        for (var c = 0; c < 3; c++) { // lol c++
          currImageData.data[idxInPixels+c] = color;
        }
      }
    }
    ctxs[1].putImageData(currImageData, 0, 0);

    enableButtons();

    var duration = +new Date() - start;
    console.log('It took '+duration+'ms to compute the FT.');
  }

  function reconstructAction() {
    // compute the h prime values
    var h_primes = [];
    var h_hats = $h();
    h_hats = Fourier.unshift(h_hats, dims);
    Fourier.invert(h_hats, h_primes);

    // store them in a nice function to match the math
    h_ = function(n, m) {
      if (arguments.length === 0) return h_primes;

      var idx = n*dims[0] + m;
      return round(h_primes[idx], 2);
    };

    // draw the pixels
    var currImageData = ctxs[2].getImageData(
      0, 0, dims[0], dims[1]
    );
    for (var n = 0; n < dims[1]; n++) {
      for (var m = 0; m < dims[0]; m++) {
        var idxInPixels = 4*(dims[0]*n + m);
        currImageData.data[idxInPixels+3] = 255; // full alpha
        for (var c = 0; c < 3; c++) { // RGB are the same, lol c++
          currImageData.data[idxInPixels+c] = h_(n, m);
        }
      }
    }
    ctxs[2].putImageData(currImageData, 0, 0);

    enableButtons();

    var duration = +new Date() - start;
    console.log(
      'It took ' + duration + 'ms to reconstruct the image.'
    );
  }

  function differenceAction() {
    // find the range of the errors
    var minError = Infinity;
    var maxError = 0;
    for (var n = 0; n < dims[1]; n++) {
      for (var m = 0; m < dims[0]; m++) {
        var error = h_(n, m) - h(n, m);
        if (error < minError) minError = error;
        if (error > maxError) maxError = error;
      }
    }

    // draw the pixels
    var currImageData = ctxs[3].getImageData(
      0, 0, dims[0], dims[1]
    );
    for (var n = 0; n < dims[1]; n++) {
      for (var m = 0; m < dims[0]; m++) {
        var idxInPixels = 4*(dims[0]*n + m);
        var error = h_(n, m) - h(n, m);
        var color = getCoolColor(
          error, [minError, maxError]
        );
        for (var c = 0; c < 3; c++) {
          currImageData.data[idxInPixels+c] = color[c];
        }
        currImageData.data[idxInPixels+3] = 255;
      }
    }
    ctxs[3].putImageData(currImageData, 0, 0);

    enableButtons();

    var duration = +new Date() - start;
    console.log(
      'It took ' + duration +
      'ms to compute the difference.'
    );
  }

  /********************
   * helper functions */
  function disableButtons(callback) {
    $s('#draw-cs-btn').disabled = true;
    $s('#draw-circle-btn').disabled = true;
    $s('#draw-grace-btn').disabled = true;
    $s('#draw-img-btn').disabled = true;
    $s('#transform-btn').disabled = true;
    $s('#reconstruct-btn').disabled = true;
    $s('#difference-btn').disabled = true;

    setTimeout(callback, 6); // 6ms for the UI to update
  }

  function enableButtons() {
    $s('#draw-cs-btn').disabled = false;
    $s('#draw-circle-btn').disabled = false;
    $s('#draw-grace-btn').disabled = false;
    $s('#draw-img-btn').disabled = false;
    $s('#transform-btn').disabled = false;
    $s('#reconstruct-btn').disabled = false;
    $s('#difference-btn').disabled = false;
  }

  // returns array of pixel colors in the image
  function getCoolColor(n, range) {
    if (n === range[0] && range[0] === range[1]) {
      return getCoolColor(2*n, [n-1, 2*n+1]);
    }

    var raw = [1.0, 1.0, 1.0]; // white

    if (n < range[0]) n = range[0];
    if (n > range[1]) n = range[1];
    var dn = range[1] - range[0];

    if (n < (range[0] + 0.25 * dn)) {
      raw[0] = 0;
      raw[1] = 4 * (n - range[0]) / dn;
    } else if (n < (range[0] + 0.5 * dn)) {
      raw[0] = 0;
      raw[2] = 1 + 4 * (range[0] + 0.25 * dn - n) / dn;
    } else if (n < (range[0] + 0.75 * dn)) {
      raw[0] = 4 * (n - range[0] - 0.5 * dn) / dn;
      raw[2] = 0;
    } else {
      raw[1] = 1 + 4 * (range[0] + 0.75 * dn - n) / dn;
      raw[2] = 0;
    }

    var color = [
      tightMap(raw[0], 0, 1, 0, 255),
      tightMap(raw[1], 0, 1, 0, 255),
      tightMap(raw[2], 0, 1, 0, 255)
    ];
    return color;
  }

  function round(n, places) {
    var mult = Math.pow(10, places);
    return Math.round(mult*n)/mult;
  }

  function tightMap(n, d1, d2, r1, r2) { // enforces boundaries
    var raw = map(n, d1, d2, r1, r2);
    if (raw < r1) return r1;
    else if (raw > r2) return r2;
    else return raw;
  }

  // given an n in [d1, d2], return a "linearly related"
  // number in [r1, r2]
  function map(n, d1, d2, r1, r2) {
    var Rd = d2-d1;
    var Rr = r2-r1;
    return (Rr/Rd)*(n - d1) + r1;
  }

  function $s(id) { // for convenience
    if (id.charAt(0) !== '#') return false;
    return document.getElementById(id.substring(1));
  }

  return {
    init: initFourierImageAnalysis
  };
})();

window.addEventListener('load', FourierImageAnalysis.init);
    
</script>