//------------------------------------------------------------------------------
#ifdef GL_ES
  precision highp float;
  precision highp int;
#endif
    #extension GL_OES_standard_derivatives : enable
//------------------------------------------------------------------------------
// Here are alll the Attributes and uniforms set by the p5.js WEBGL mode architecture
// https://github.com/processing/p5.js/blob/master/developer_docs/webgl_mode_architecture.md
//------------------------------------------------------------------------------
// attributes: these come from the gl context memoryBarrierBuffer
// Names are defined by p5 webgl implementation
//Geometry
attribute vec3 aPosition;
attribute vec3 aNormal;
attribute vec4 aDirection;
attribute vec2 aTexCoord;

// color
attribute vec4 aVertexColor;

//------------------------------------------------------------------------------
// Uniforms are variable passed from the program to the shader. In this case, javascript
uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;
uniform vec4 uViewPort;

// geometry
uniform mat3 uNormalMatrix;
uniform float uStrokeWeight;

// Color
uniform vec4 uMaterialColor;

// Light Parameters
uniform int uAmbientLightCount;
uniform int uDirectionalLightCount;
uniform int uPointLightCount;
uniform vec3 uAmbientColor[8];
uniform vec3 uLightingDirection[8];
uniform vec3 uDirectionalColor[8];
uniform vec3 uPointLightLocation[8];
uniform vec3 uPointLightColor[8];
uniform bool uSpecular;
uniform int uShininess;
uniform bool uUseLighting;

//Texture
uniform sampler2D uSampler;
uniform bool isTexture;

// general
uniform float uResolution;
uniform float uPointSize;

// scaling
uniform float zscale;
//------------------------------------------------------------------------------
// varying: these go to the frag shader out
varying vec3 var_vertPos;
varying vec4 var_vertCol;
varying vec3 var_vertNormal;
varying vec2 var_vertTexCoord;
varying vec4 v_color;
//------------------------------------------------------------------------------
void main()
{
  vec4 temp = vec4(aPosition, 1.0);
  temp.z *= zscale;
  gl_Position = uProjectionMatrix * uModelViewMatrix * temp;
  // gl_Position = vec4(aPosition, 1.0);
  v_color = gl_Position * 0.5 + 0.5;
  // Pass through to frag shader
  var_vertPos      = aPosition;
  var_vertCol      = aVertexColor;
  var_vertNormal   = aNormal;
  var_vertTexCoord = aTexCoord;
}
