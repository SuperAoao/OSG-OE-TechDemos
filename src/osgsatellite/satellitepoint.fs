varying int window_pos_x;
varying int window_pos_y;

varying float n_pe;
varying int ptSize;
varying int highlightType;
uniform float n_e;  // 每1个灰度值对应的电子数
void main()
{
    gl_FragColor = vec4(gl_FragCoord.xyz, 1.0);//vec4(1.0, 0.0, 0.0, 1.0); // Set the color of the point (red)

    if (gl_FragCoord.x <= 1.0)
    {
        gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
    }
}
