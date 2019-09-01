package ioio.examples.hello;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Timer;
import java.util.TimerTask;

//import com.thousandthoughts.tutorials.SensorFusionActivity.calculateFusedOrientationTask;
import ioio.lib.util.BaseIOIOLooper;
import ioio.lib.util.IOIOLooper;
import ioio.lib.util.android.IOIOActivity;
import ioio.lib.api.DigitalInput;
import ioio.lib.api.DigitalOutput;
import ioio.lib.api.PulseInput;
import ioio.lib.api.PwmOutput;
import ioio.lib.api.exception.ConnectionLostException;
import android.content.Context;
import android.hardware.Sensor;
import android.hardware.SensorEvent;
import android.hardware.SensorEventListener;
import android.hardware.SensorManager;
import android.media.MediaScannerConnection;
import android.os.Build;
import android.os.Bundle;
import android.util.Log;
import android.widget.TextView;
import android.widget.ToggleButton;
import android.annotation.SuppressLint;
import android.annotation.TargetApi;
import android.app.Activity;

/**
 * This is the main activity of the HelloIOIO example application.
 * 
 * It displays a toggle button on the screen, which enables control of the
 * on-board LED. This example shows a very simple usage of the IOIO, by using
 * the {@link IOIOActivity} class. For a more advanced use case, see the
 * HelloIOIOPower example.
 */

public class MainActivity extends IOIOActivity implements SensorEventListener{
	
	int igravado=0;
	double Ttotal=0;
	int ngravador=2000;
	
	float xl=0;
	float yl=0;
	float zl=0;
	
	
	public double[] MemoY=new double[100];
	int py=0;
	int tamanhobuff=10;
	public static double saidamedia=0;
	
	float[] gravity = new float[3];
	
	float[] Rx = new float[9];
	float[] Ry = new float[9];
	float[] Rz = new float[9];
	float[] Rt = new float[3];
	
	float[] RxRy= new float[9];
	float[] RxRyRz= new float[9];
	

	public long time=0;
	public boolean timeiniciador=false;
	public long timem1=0;
	
	public float sinalnormalizado;
	public float sinal;
	public float sinalpulso;
	public float sinalpulso2;
	
	public float sinalseno=0;
	public float amplitudeseno=0;
	public float offset1=0;
	public boolean sentido1=false;
	public float sinalseno2=0;
	public float amplitudeseno2=0;
	public float offset2=0;
	public boolean sentido2=false;
	public int ktraj=1;
	
	
	
	
	
	public static float Ak=0;
    public static float Vk=0;
    public static float Vkm1=0;
    public float Xk=0;
    public float Akm1;
    public float Xkm1=0;
    public float Dt=0;

    public static float dt;
    
	
	
	private ToggleButton button_;
	private SensorManager mSensorManager = null;
	
	private float[] Xkgravado = new float[ngravador+1];
	private float[] Ykgravado = new float[ngravador+1];
	private float[] Zkgravado = new float[ngravador+1];
	
	private float[] OriZ= new float[ngravador+1];
	private float[] OriY= new float[ngravador+1];
	private float[] OriX= new float[ngravador+1];
	private double[] DT= new double[ngravador+1];
	
	
	
	//Velocidades angulares do giro
    private float[] gyro = new float[3];
 
    // Matriz de rota��o dos dados do giro
    private float[] gyroMatrix = new float[9];
 
    // Angulos de orienta��o da matriz do giroo
    private float[] gyroOrientation = new float[3];
 
    // Vetor de campo magn�tico
    private float[] magnet = new float[3];
 
    // Vetor de Acelera��o
    private float[] accel = new float[3];
 
    // orientation angles from accel and magnet
    private float[] accMagOrientation = new float[3];
 
    // final orientation angles from sensor fusion
    private float[] fusedOrientation = new float[3];
 
    // accelerometer and magnetometer based rotation matrix
    private float[] rotationMatrix = new float[9];
    
    public static final float EPSILON = 0.000000001f;
    private static final float NS2S = 1.0f / 1000000000.0f;
	private long timestamp;
	private boolean initState = true;
	public static final int TIME_CONSTANT = 20;
	public static final float FILTER_COEFFICIENT = 0.90f;
	private Timer fuseTimer = new Timer();
	
	
	
	
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.main);
		button_ = (ToggleButton) findViewById(R.id.button);
		mSensorManager = (SensorManager) this.getSystemService(SENSOR_SERVICE);
		//fuseTimer.scheduleAtFixedRate(new calculateFusedOrientationTask(),
          //      1000, TIME_CONSTANT);
		
		gyroOrientation[0] = 0.0f;
        gyroOrientation[1] = 0.0f;
        gyroOrientation[2] = 0.0f;
 
        // initialise gyroMatrix with identity matrix
        gyroMatrix[0] = 1.0f; gyroMatrix[1] = 0.0f; gyroMatrix[2] = 0.0f;
        gyroMatrix[3] = 0.0f; gyroMatrix[4] = 1.0f; gyroMatrix[5] = 0.0f;
        gyroMatrix[6] = 0.0f; gyroMatrix[7] = 0.0f; gyroMatrix[8] = 1.0f;
		
		
        initListeners();
		
        fuseTimer.scheduleAtFixedRate(new calculateFusedOrientationTask(),
                500, TIME_CONSTANT);
		
	}

	@Override
    public void onStop() {
    	super.onStop();
    	// unregister sensor listeners to prevent the activity from draining the device's battery.
    	mSensorManager.unregisterListener(this);
    }
	
    @Override
    protected void onPause() {
        super.onPause();
        // unregister sensor listeners to prevent the activity from draining the device's battery.
        mSensorManager.unregisterListener(this);
    }
    
    @Override
    public void onResume() {
    	super.onResume();
    	// restore the sensor listeners when user resumes the application.
    	initListeners();
    }
    
    public void initListeners(){
        mSensorManager.registerListener(this,
            mSensorManager.getDefaultSensor(Sensor.TYPE_ACCELEROMETER),
            SensorManager.SENSOR_DELAY_FASTEST);
     
        mSensorManager.registerListener(this,
            mSensorManager.getDefaultSensor(Sensor.TYPE_GYROSCOPE),
            SensorManager.SENSOR_DELAY_FASTEST);
     
        mSensorManager.registerListener(this,
            mSensorManager.getDefaultSensor(Sensor.TYPE_MAGNETIC_FIELD),
            SensorManager.SENSOR_DELAY_FASTEST);
    }
	
	
	
	
	@Override
	protected IOIOLooper createIOIOLooper() {
		return new Looper();
	}

	@Override
	public void onAccuracyChanged(Sensor arg0, int arg1) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void onSensorChanged(SensorEvent event) {
		
		
		
switch(event.sensor.getType()) {
	    
		
		case Sensor.TYPE_ACCELEROMETER:
	        
			// Copia os dados do aceler�metro para um vetor accel e calcula orienta��o
	        System.arraycopy(event.values, 0, accel, 0, 3);
	        
			
	        time=event.timestamp;
			
			if(!timeiniciador)
			{
				timem1=time;
				timeiniciador=true;
			}
	        
	        
	        calculateAccMagOrientation();
	        estimardeslocamento();
	        
	        break;
	 
	    case Sensor.TYPE_GYROSCOPE:
	        // Processa os dados do giroscopio
	        gyroFunction(event);
	        break;
	 
	    case Sensor.TYPE_MAGNETIC_FIELD:
	        // Copia dados do magnetometro para um vetor
	        System.arraycopy(event.values, 0, magnet, 0, 3);
	        break;
	    }
		
	}
	
	
	private void getRotationVectorFromGyro(float[] gyroValues,
            float[] deltaRotationVector,
            float timeFactor)
	{
		float[] normValues = new float[3];
		
		// Calculate the angular speed of the sample
		float omegaMagnitude =
		(float)Math.sqrt(gyroValues[0] * gyroValues[0] +
		gyroValues[1] * gyroValues[1] +
		gyroValues[2] * gyroValues[2]);
		
		// Normalize the rotation vector if it's big enough to get the axis
		if(omegaMagnitude > EPSILON) {
		normValues[0] = gyroValues[0] / omegaMagnitude;
		normValues[1] = gyroValues[1] / omegaMagnitude;
		normValues[2] = gyroValues[2] / omegaMagnitude;
		}
		
		// Integrate around this axis with the angular speed by the timestep
		// in order to get a delta rotation from this sample over the timestep
		// We will convert this axis-angle representation of the delta rotation
		// into a quaternion before turning it into the rotation matrix.
		float thetaOverTwo = omegaMagnitude * timeFactor;
		float sinThetaOverTwo = (float)Math.sin(thetaOverTwo);
		float cosThetaOverTwo = (float)Math.cos(thetaOverTwo);
		deltaRotationVector[0] = sinThetaOverTwo * normValues[0];
		deltaRotationVector[1] = sinThetaOverTwo * normValues[1];
		deltaRotationVector[2] = sinThetaOverTwo * normValues[2];
		deltaRotationVector[3] = cosThetaOverTwo;
	}
	private float[] matrixMultiplication(float[] A, float[] B) {
        float[] result = new float[9];
     
        result[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
        result[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
        result[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
     
        result[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
        result[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
        result[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
     
        result[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
        result[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
        result[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
     
        return result;
    }
	private float[] getRotationMatrixFromOrientation(float[] o) {
        float[] xM = new float[9];
        float[] yM = new float[9];
        float[] zM = new float[9];
     
        float sinX = (float)Math.sin(o[1]);
        float cosX = (float)Math.cos(o[1]);
        float sinY = (float)Math.sin(o[2]);
        float cosY = (float)Math.cos(o[2]);
        float sinZ = (float)Math.sin(o[0]);
        float cosZ = (float)Math.cos(o[0]);
       
        // Rota��o em torno de X (pitch)
        xM[0] = 1.0f; xM[1] = 0.0f; xM[2] = 0.0f;
        xM[3] = 0.0f; xM[4] = cosX; xM[5] = sinX;
        xM[6] = 0.0f; xM[7] = -sinX; xM[8] = cosX;
     
        // Rota��o em torno de Y (roll)
        yM[0] = cosY; yM[1] = 0.0f; yM[2] = sinY;
        yM[3] = 0.0f; yM[4] = 1.0f; yM[5] = 0.0f;
        yM[6] = -sinY; yM[7] = 0.0f; yM[8] = cosY;
     
        // Rota��o em torno de Z (azimuth)
        zM[0] = cosZ; zM[1] = sinZ; zM[2] = 0.0f;
        zM[3] = -sinZ; zM[4] = cosZ; zM[5] = 0.0f;
        zM[6] = 0.0f; zM[7] = 0.0f; zM[8] = 1.0f;
     
        // Ordem de rota��o y, x, z (roll, pitch, azimuth)
        float[] resultMatrix = matrixMultiplication(xM, yM);
        resultMatrix = matrixMultiplication(zM, resultMatrix);
        return resultMatrix;
    }
	public void gyroFunction(SensorEvent event) {
        // don't start until first accelerometer/magnetometer orientation has been acquired
        if (accMagOrientation == null)
            return;
     
        // initialisation of the gyroscope based rotation matrix
        if(initState) {
            float[] initMatrix = new float[9];
            initMatrix = getRotationMatrixFromOrientation(accMagOrientation);
            float[] test = new float[3];
            SensorManager.getOrientation(initMatrix, test);
            gyroMatrix = matrixMultiplication(gyroMatrix, initMatrix);
            initState = false;
        }
     
        // copy the new gyro values into the gyro array
        // convert the raw gyro data into a rotation vector
        float[] deltaVector = new float[4];
        if(timestamp != 0) {
            final float dT = (event.timestamp - timestamp) * NS2S;
        System.arraycopy(event.values, 0, gyro, 0, 3);
        getRotationVectorFromGyro(gyro, deltaVector, dT / 2.0f);
        }
     
        // measurement done, save current time for next interval
        timestamp = event.timestamp;
     
        // convert rotation vector into rotation matrix
        float[] deltaMatrix = new float[9];
        SensorManager.getRotationMatrixFromVector(deltaMatrix, deltaVector);
     
        // apply the new rotation interval on the gyroscope based rotation matrix
        gyroMatrix = matrixMultiplication(gyroMatrix, deltaMatrix);
     
        // get the gyroscope based orientation from the rotation matrix
        SensorManager.getOrientation(gyroMatrix, gyroOrientation);
        
        
        
        
        
        
    }
	 

	class calculateFusedOrientationTask extends TimerTask {
        public void run() {
            float oneMinusCoeff = 1.0f - FILTER_COEFFICIENT;
            
            /*
             * Fix for 179� <--> -179� transition problem:
             * Check whether one of the two orientation angles (gyro or accMag) is negative while the other one is positive.
             * If so, add 360� (2 * math.PI) to the negative value, perform the sensor fusion, and remove the 360� from the result
             * if it is greater than 180�. This stabilizes the output in positive-to-negative-transition cases.
             */
            
            // azimuth
            if (gyroOrientation[0] < -0.5 * Math.PI && accMagOrientation[0] > 0.0) {
            	fusedOrientation[0] = (float) (FILTER_COEFFICIENT * (gyroOrientation[0] + 2.0 * Math.PI) + oneMinusCoeff * accMagOrientation[0]);
        		fusedOrientation[0] -= (fusedOrientation[0] > Math.PI) ? 2.0 * Math.PI : 0;
            }
            else if (accMagOrientation[0] < -0.5 * Math.PI && gyroOrientation[0] > 0.0) {
            	fusedOrientation[0] = (float) (FILTER_COEFFICIENT * gyroOrientation[0] + oneMinusCoeff * (accMagOrientation[0] + 2.0 * Math.PI));
            	fusedOrientation[0] -= (fusedOrientation[0] > Math.PI)? 2.0 * Math.PI : 0;
            }
            else {
            	fusedOrientation[0] = FILTER_COEFFICIENT * gyroOrientation[0] + oneMinusCoeff * accMagOrientation[0];
            }
            
            // pitch
            if (gyroOrientation[1] < -0.5 * Math.PI && accMagOrientation[1] > 0.0) {
            	fusedOrientation[1] = (float) (FILTER_COEFFICIENT * (gyroOrientation[1] + 2.0 * Math.PI) + oneMinusCoeff * accMagOrientation[1]);
        		fusedOrientation[1] -= (fusedOrientation[1] > Math.PI) ? 2.0 * Math.PI : 0;
            }
            else if (accMagOrientation[1] < -0.5 * Math.PI && gyroOrientation[1] > 0.0) {
            	fusedOrientation[1] = (float) (FILTER_COEFFICIENT * gyroOrientation[1] + oneMinusCoeff * (accMagOrientation[1] + 2.0 * Math.PI));
            	fusedOrientation[1] -= (fusedOrientation[1] > Math.PI)? 2.0 * Math.PI : 0;
            }
            else {
            	fusedOrientation[1] = FILTER_COEFFICIENT * gyroOrientation[1] + oneMinusCoeff * accMagOrientation[1];
            }
            
            // roll
            if (gyroOrientation[2] < -0.5 * Math.PI && accMagOrientation[2] > 0.0) {
            	fusedOrientation[2] = (float) (FILTER_COEFFICIENT * (gyroOrientation[2] + 2.0 * Math.PI) + oneMinusCoeff * accMagOrientation[2]);
        		fusedOrientation[2] -= (fusedOrientation[2] > Math.PI) ? 2.0 * Math.PI : 0;
            }
            else if (accMagOrientation[2] < -0.5 * Math.PI && gyroOrientation[2] > 0.0) {
            	fusedOrientation[2] = (float) (FILTER_COEFFICIENT * gyroOrientation[2] + oneMinusCoeff * (accMagOrientation[2] + 2.0 * Math.PI));
            	fusedOrientation[2] -= (fusedOrientation[2] > Math.PI)? 2.0 * Math.PI : 0;
            }
            else {
            	fusedOrientation[2] = FILTER_COEFFICIENT * gyroOrientation[2] + oneMinusCoeff * accMagOrientation[2];
            }
     
           
            //gyroMatrix = getRotationMatrixFromOrientation(fusedOrientation);
            //System.arraycopy(fusedOrientation, 0, gyroOrientation, 0, 3);
            
        }
    }
	
	        
	public void gravardados(){
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "Ay.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(Ykgravado[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
		
		catch (IOException e)
		    {
		    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
		    }
		
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "Ax.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(Xkgravado[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "Az.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(Zkgravado[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "OriZ.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(OriZ[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "OriY.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(OriY[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "OriX.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(OriX[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		try
		   {
			int i=0;
		                            
		   File traceFile = new File(this.getExternalFilesDir(null), "DT.txt");
		   if (!traceFile.exists())
		      traceFile.createNewFile();
		                        
		   BufferedWriter writer = new BufferedWriter(new FileWriter(traceFile, true /*append*/));
		   
		   while(i<ngravador)
		   {
		   writer.write("; "+Double.toString(DT[i]));
		   i++;
		   }
		   
		   i=0;
		   writer.close();
		                          
		    MediaScannerConnection.scanFile(this,
		                                     new String[] { traceFile.toString() },
		                                     null,
		                                     null);
		    }
	
		catch (IOException e)
	    {
	    //Log.e("com.cindypotvin.FileTest", "Unable to write to the TraceFile.txt file.");
	    }
		
		igravado=0;
		
		
		
		}
		

	
	 
	public void estimardeslocamento(){
		float alpha=(float) 0.9999;
		gravity[0] = alpha * gravity[0] + (1 - alpha) * accel[0];
		  gravity[1] = alpha * gravity[1] + (1 - alpha) * accel[1];
		  gravity[2] = alpha * gravity[2] + (1 - alpha) * accel[2];

		  // Remove the gravity contribution with the high-pass filter.
		  xl = accel[0]; 
		  yl = accel[1]; 
		  zl = accel[2];
		  
		TextView desloc=(TextView)findViewById(R.id.desloc);
		TextView vari=(TextView)findViewById(R.id.var);
		
		
		desloc.setText("Deslocamento "+Double.toString(Xk*10));
		
		
		Xkgravado[igravado]= xl;
		Ykgravado[igravado]= yl;
		Zkgravado[igravado]= zl;
		
		if(gyroOrientation[0]<0)
        {
        	gyroOrientation[0]=(float) (gyroOrientation[0]+2*(Math.PI));
        }
		
		OriZ[igravado]=(float)(gyroOrientation[0]*180/Math.PI);
		OriY[igravado]=(float)(gyroOrientation[1]*180/Math.PI);
		OriX[igravado]=(float)(gyroOrientation[2]*180/Math.PI);
		DT[igravado]=(float)(Ttotal);
		
		
		//Xkgravado[igravado]=(float) (Xk*10);
		if(igravado==ngravador)
		{
		gravardados();	
		}
		igravado++;
		
		}
		
	
	        
	public void calculateAccMagOrientation() {
		
		TextView acelx=(TextView)findViewById(R.id.acelx);
		TextView acely=(TextView)findViewById(R.id.acely);
		TextView acelz=(TextView)findViewById(R.id.acelz);
		
		TextView Rtx=(TextView)findViewById(R.id.Rtx);
		TextView Pitch=(TextView)findViewById(R.id.Y);
		TextView Yaw=(TextView)findViewById(R.id.Z);
		TextView Roll=(TextView)findViewById(R.id.Roll);
		
	    if(SensorManager.getRotationMatrix(rotationMatrix, null, accel, magnet)) {
	        SensorManager.getOrientation(rotationMatrix, accMagOrientation);
	        
	       
	        Rz[0]=(float) Math.cos(accMagOrientation[0]) ; Rz[1]=(float) Math.sin(accMagOrientation[0]);Rz[2]= 0;
	        Rz[3]=(float) -Math.sin(accMagOrientation[0]);Rz[4]= (float) Math.cos(accMagOrientation[0]);Rz[5]= 0;
	        Rz[6]=0; Rz[7]=0;Rz[8]=1;
			
			Ry[0]=(float) Math.cos(accMagOrientation[1]); Ry[1]=0; Ry[2]=(float) -Math.sin(accMagOrientation[1]);
	        Ry[3]= 0; Ry[4]= 1; Ry[5]=0;
	        Ry[6]=(float) Math.sin(accMagOrientation[1]);Ry[7]=0; Ry[8]=(float) Math.cos(accMagOrientation[1]);
			
	        Rx[0]=1; Rx[1]=0; Rx[2]=0;
	        Rx[3]= 0; Rx[4]= (float) Math.cos(accMagOrientation[2]); Rx[5]=(float) Math.sin(accMagOrientation[2]);
	        Rx[6]=0;Rx[7]=(float) -Math.sin(accMagOrientation[2]); Rx[8]=(float) Math.cos(accMagOrientation[2]);
			
			
	        
	        //accMagOrientation[0]=(float) ((accMagOrientation[0]*180)/Math.PI);
	        //accMagOrientation[1]=(float) ((accMagOrientation[1]*180)/Math.PI);
	        //accMagOrientation[2]=(float) ((accMagOrientation[2]*180)/Math.PI);
	        
	        /*
	        Rz[0]=(float) Math.cos(fusedOrientation[0]) ; Rz[1]=(float) Math.sin(fusedOrientation[0]);Rz[2]= 0;
	        Rz[3]=(float) -Math.sin(fusedOrientation[0]);Rz[4]= (float) Math.cos(fusedOrientation[0]);Rz[5]= 0;
	        Rz[6]=0; Rz[7]=0;Rz[8]=1;
			
			Ry[0]=(float) Math.cos(fusedOrientation[1]); Ry[1]=0; Ry[2]=(float) -Math.sin(fusedOrientation[1]);
	        Ry[3]= 0; Ry[4]= 1; Ry[5]=0;
	        Ry[6]=(float) Math.sin(fusedOrientation[1]);Ry[7]=0; Ry[8]=(float) Math.cos(fusedOrientation[1]);
			
	        Rx[0]=1; Rx[1]=0; Rx[2]=0;
	        Rx[3]= 0; Rx[4]= (float) Math.cos(fusedOrientation[2]); Rx[5]=(float) Math.sin(fusedOrientation[2]);
	        Rx[6]=0;Rx[7]=(float) -Math.sin(fusedOrientation[2]); Rx[8]=(float) Math.cos(fusedOrientation[2]);
	        */
	        
	        RxRy = matrixMultiplication(Rx, Ry);
			RxRyRz= matrixMultiplication(RxRy,Rz);
	        
			Rt[0]= (float) (RxRyRz[2]*9.81);
			Rt[1]= (float)(RxRyRz[5]*9.81);
			Rt[2]= (float) (RxRyRz[8]*10.6);
	        
	        
			acelx.setText("Ax "+Float.toString(accel[0]+Rt[1]));
	        acely.setText("Ay "+Float.toString(accel[1]));
	        acelz.setText("Az "+Float.toString(accel[2]));
	        
	        
	        Dt= (time-timem1)*NS2S;
	        timem1=time;
	        Ttotal=Dt+Ttotal;
	       
	        // Botar fusedorientation aqui
	        Rtx.setText("T amostragem "+Float.toString(Dt));       
	        Yaw.setText(Float.toString((float)   (fusedOrientation[0]*180/Math.PI)));  // Z
	        Pitch.setText(Float.toString((float) (fusedOrientation[1]*180/Math.PI))); // Y
	        Roll.setText(Float.toString((float)  (fusedOrientation[2]*180/Math.PI))); //  X
	        
	    }
	}
	
	
	
	
	
	
	
	/**
	 * This is the thread on which all the IOIO activity happens. It will be run
	 * every time the application is resumed and aborted when it is paused. The
	 * method setup() will be called right after a connection with the IOIO has
	 * been established (which might happen several times!). Then, loop() will
	 * be called repetitively until the IOIO gets disconnected.
	 */
	class Looper extends BaseIOIOLooper {

		private DigitalOutput led_;
		private DigitalOutput saidaseg;
		private DigitalOutput saidadig;
		private DigitalOutput saidadig2;
		private PwmOutput pwm;
		private PwmOutput pwm2;
		private PulseInput pulso;
		private PulseInput pulso2;
		public float k=0;
		public float k2=0;
		
		
		/**
		 * Called every time a connection with IOIO has been established.
		 * Typically used to open pins.
		 * 
		 * @throws ConnectionLostException
		 *             When IOIO connection is lost.
		 * 
		 * @see ioio.lib.util.AbstractIOIOActivity.IOIOThread#setup()
		 */
		@Override
		protected void setup() throws ConnectionLostException {
			led_ = ioio_.openDigitalOutput(0, true);
			saidaseg= ioio_.openDigitalOutput(1,true);
			saidadig= ioio_.openDigitalOutput(3,true);
			saidadig2= ioio_.openDigitalOutput(4,true);
			pwm = ioio_.openPwmOutput(36, 100);
			pwm2= ioio_.openPwmOutput(37, 100);
			pulso = ioio_.openPulseInput(new DigitalInput.Spec(7),
	                 PulseInput.ClockRate.RATE_250KHz,
	                 PulseInput.PulseMode.POSITIVE,
	                  false);	
			
			pulso2 = ioio_.openPulseInput(new DigitalInput.Spec(6),
	                 PulseInput.ClockRate.RATE_250KHz,
	                 PulseInput.PulseMode.POSITIVE,
	                  false);	
			
		}

		/**
		 * Called repetitively while the IOIO is connected.
		 * 
		 * @throws ConnectionLostException
		 *             When IOIO connection is lost.
		 * 
		 * @see ioio.lib.util.AbstractIOIOActivity.IOIOThread#loop()
		 */
		@Override
		public void loop() throws ConnectionLostException {
			
			//saidaseg.write(true);
			saidadig.write(false);
			saidadig.write(false);
			pwm.setDutyCycle((float) (1)); //1 � o m�nimo, 0 � o m�ximo
			pwm2.setDutyCycle((float) (1));
			
			led_.write(!button_.isChecked());
			
			
			while(button_.isChecked())
			{

			//Entradas


			/*if(ktraj==1)                    // Modificação 13/07/17 - 18hr
			{
				amplitudeseno=(float) 0.9;    //0.5  reta - pwm2(cadeira de rodas)
				amplitudeseno2=(float) 0.87;  //0.43 reta - pwm1(cadeira de rodas)
			}*/

////////////////////////////////////////////////----------------------
			/*	if(ktraj==1)                    // Modificação 13/07/17 - 15hr
				{
					amplitudeseno=(float) 0.52;   //0.5 reta    pwm2(cadeira de rodas)
					amplitudeseno2=(float) 0.4;  //0.43 reta     pwm1(cadeira de rodas)
				}*/
/////////////////////////////////////////////////----------------------
			/*if(ktraj==1)
			{
			amplitudeseno=(float) 0.5;   //0.5 reta
			amplitudeseno2=(float) 0.6;  //0.43 reta
			}
////////////////////////////////////////////////-----------------------
			if(ktraj==2)
			{
				amplitudeseno=(float) 0.8;   //0.5 reta
				amplitudeseno2=(float) 0.0;  //0.43 reta
			}
////////////////////////////////////////////////-----------------------
			if(ktraj==3)
			{
			amplitudeseno=(float) 0.4;   //0.5 reta
			amplitudeseno2=(float) 0.6;  //0.43 reta
			}

////////////////////////////////////////////////-----------------------
			if(ktraj==4)
			{
				amplitudeseno=(float) 0;   //0.5 reta
				amplitudeseno2=(float) 0.73;  //0.43 reta
				ktraj=0;
			}*/
			
			
			
			//Se sentido=true, offsetdesligado
			sentido1=false;  
			sentido2=false;
			
			
			
			
			//Ajuste de sinal
			saidadig.write(sentido1); //true = Offset desativado, false = offset ativado
			saidadig2.write(sentido2);
			
			
			if(sentido1==false)
			{
			offset1=(1-amplitudeseno)*2;
			sinalseno=(float) ((amplitudeseno)*(Math.sin(2*k*Math.PI/180)+Math.sin(90*Math.PI/180))+offset1);
			sinalseno=sinalseno/2;
			}
			else
			{
			sinalseno=(float) ((amplitudeseno)*(Math.sin(2*k*Math.PI/180)+Math.sin(90*Math.PI/180))+offset1);
				sinalseno=sinalseno/2;
			}
			
			
			if(sentido2==false)
			{
			offset2=(1-amplitudeseno2)*2;
			sinalseno2=(float) ((amplitudeseno2)*(Math.sin(2*k*Math.PI/180)+Math.sin(90*Math.PI/180))+offset2);
			sinalseno2=sinalseno2/2;
			}
			else
			{
			sinalseno2=(float) ((amplitudeseno2)*(Math.sin(2*k*Math.PI/180)+Math.sin(90*Math.PI/180))+offset2);
				sinalseno2=sinalseno2/2;
			}
			
			
			
			
			//Trejet�ria
			pwm.setDutyCycle(sinalseno);  //0 a 2V 2.5 a 4.5
			pwm2.setDutyCycle(sinalseno2);
			
			
			if(ktraj>269)
			{
			ktraj=ktraj+1;
			}
			
			if(k>359)
			{
				
				k=0;
			}
			
			k++;
			
			
			
			
			
			
			
			
			
			
			
			

			try {
				Thread.sleep(10);
			} catch (InterruptedException e) {
			}
			
			
		}
			
			
			
		}
	}
}

	/**
	 * A method to create our IOIO thread.
	 * 
	 * @see ioio.lib.util.AbstractIOIOActivity#createIOIOThread()
	 */
	
