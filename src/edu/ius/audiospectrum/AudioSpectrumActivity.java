/*
 *     Copyright (C) 2013  Raymond Wisman
 * 			Indiana University SE
 * 			December 29, 2013
 * 
  	AudioSpectrum displays the frequency spectrum of a WAV file or audio sourced from AudioTime or AudioTime+.

	The spectrum is displayed in linear and logarithmic scales.
	Panning and zooming allow exploration of the fine details of the spectrum.
	Ensemble spectral averaging is used to construct the spectrum of continuous audio too large for a single FFT transformation.

	The application is designed for use in science education experiments that:

		investigate the frequency components of sound,
		determine the dominant frequency of audio.
		
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

package edu.ius.audiospectrum;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.text.DecimalFormat;
import java.text.FieldPosition;
import java.text.Format;
import java.text.ParsePosition;

import android.app.Activity;
import android.app.AlertDialog;
import android.app.Dialog;
import android.app.ProgressDialog;
import android.content.ContentResolver;
import android.content.DialogInterface;
import android.content.Intent;
import android.graphics.Color;
import android.graphics.DashPathEffect;
import android.graphics.PointF;
import android.net.Uri;
import android.os.Bundle;
import android.os.Handler;
import android.os.Message;
import android.util.Pair;
import android.view.GestureDetector;
import android.view.GestureDetector.OnDoubleTapListener;
import android.view.MotionEvent;
import android.view.View;
import android.view.GestureDetector.OnGestureListener;
import android.widget.AdapterView;
import android.widget.ArrayAdapter;
import android.widget.ImageButton;
import android.widget.ListView;

import com.androidplot.LineRegion;
import com.androidplot.ui.YLayoutStyle;
import com.androidplot.ui.YPositionMetric;
import com.androidplot.util.PixelUtils;
import com.androidplot.xy.*;

import edu.ius.fftresult.FFTresult;
import edu.ius.fft.FFT;

public class AudioSpectrumActivity extends Activity implements OnGestureListener, OnDoubleTapListener {
	private static final int SPECTRUMSIZE = 400;
	private static final int MAXSAMPLES = 65536;
	private static final long ELAPSEDTIME = 50l;
	private static final int ABOUT = 0;
	private static final int INFO = 1;

	private final Activity activity = this;
	private static final String[] DIALOGMENU = new String[] { "About", "Info" };

	private static XYPlot plot;
	private static Pair<Integer, XYSeries> selection;
	private static XValueMarker xMarker = null;
	private static SimpleXYSeries seriesData = null;
	private static BarFormatter seriesFormat = null;
	private static GraphXLabelFormat graphXLabelFormat = null;

	private static double magnitudeFFT[] = null;
	private static int frequencyFFT[]=null;
	private static int sampleRate;
	private static int samples;
	private static int samplesToDisplay;
	private static float maxF, minF, maxX;
	private static boolean isReloading = false;
	private static boolean isDrawing = false;
	private static boolean scaleLinear;
	
	private static long lastTime=System.currentTimeMillis();
	private static Thread threadReload=null, threadRedraw=null;
	private static GestureDetector detector;
	
	private ProgressThread progThread;
    private ProgressDialog progDialog;
    
    private Uri uri = null;

	@SuppressWarnings("deprecation")
	@Override
	public void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		
		scaleLinear = true;
		
		setContentView(R.layout.audiospectrum);

		detector = new GestureDetector(this, this);

        seriesFormat = new BarFormatter(Color.LTGRAY, Color.LTGRAY);

		seriesData = new SimpleXYSeries("Hz");

		plot = (XYPlot) findViewById(R.id.spectrumPlot);
		plot.getGraphWidget().setRangeValueFormat(new DecimalFormat("#.##"));
		plot.setTicksPerRangeLabel(1);
		plot.setRangeLowerBoundary(0, BoundaryMode.FIXED);
		plot.getGraphWidget().setRangeLabelOrientation(-45);

		plot.setTicksPerDomainLabel(1);
		plot.getGraphWidget().setDomainLabelOrientation(-45);
		
		graphXLabelFormat = new GraphXLabelFormat(scaleLinear);
		plot.getGraphWidget().setDomainValueFormat(graphXLabelFormat);
		
		startReloadThread();
		startRedrawThread();
		setButtonHandlers();
		
 		Intent intent = getIntent(); 
 		String action = intent.getAction(); 
 
//		loadInitialSeries();
 		
		if(action.equals("edu.ius.audiospectrum.ACTION") &&  intent.getExtras() != null) {			// Started from another app passing FFTresult object file name
    		try {
    			String filename = intent.getExtras().getString("FFTresult");
    			ObjectInputStream in = new ObjectInputStream(new FileInputStream(filename));
       			FFTresult fftResult = (FFTresult)in.readObject();
//    			Log.d("3 ****edu.ius.AudioSpectrum.ACTION onCreate ",filename+" " + fftResult);
       			graphSpectrum(fftResult.magnitudeFFT, fftResult.frequencyFFT, fftResult.sampleRate, fftResult.samples);
				in.close();
				(new File(filename)).delete();
	   		} catch (ClassNotFoundException e) {
    			e.printStackTrace();
	   		} catch (FileNotFoundException e) {
    			e.printStackTrace();
    		} catch (IOException e) {
    			e.printStackTrace();
    		}
		}
		else {
			intent = getIntent(); 							
			action = intent.getAction(); 
			if(!action.equals(Intent.ACTION_VIEW))  {			// Started directly as app
			 	intent = new Intent();
		        intent.addCategory(Intent.CATEGORY_OPENABLE);
 		        intent.setType("audio/wav");
		        intent.setAction(Intent.ACTION_GET_CONTENT);
	        	startActivityForResult(Intent.createChooser(intent, "Select a WAV file for "+getString(R.string.app_name)+" via a file manager."),1);	
			}
			else {
				uri = getIntent().getData();
				showDialog(0);									// Start ProgressBar, input WAV file and run FFT
			}
				
		}
	}
	
	@Override
	protected void onActivityResult(int requestCode, int resultCode, Intent data) {
		if (resultCode == RESULT_OK && requestCode == 1) {
			uri = data.getData();
			showDialog(0);
		}
		else
			errorDialog(getString(R.string.app_name)+" requires a file manager to select a WAV file or AudioTime or AudioTime+ apps (View the FFT spectrum).");
	}

	@Override
    protected Dialog onCreateDialog(int id) {
        switch(id) {
        case 0:                      // Spinner
            progDialog = new ProgressDialog(this);
            progDialog.setCanceledOnTouchOutside(false);
            progDialog.setProgressStyle(ProgressDialog.STYLE_SPINNER);
            progDialog.setMessage("Computing FFT...");
            progThread = new ProgressThread(handler);
            progThread.start();
            return progDialog;
        default:
            return null;
        }
    }
    
    final Handler handler = new Handler() {
        @SuppressWarnings("deprecation")
		public void handleMessage(Message msg) {
            dismissDialog(0);
        }
    };
 	
    private class ProgressThread extends Thread {	
        Handler mHandler;
    
        ProgressThread(Handler h) {
            mHandler = h;
        }
        
        @Override
        public void run() {   
        	loadInitialSeries();
            Message msg = mHandler.obtainMessage();
            mHandler.sendMessage(msg);
        }
     }
    
    private void errorDialog(final String msg){
	    runOnUiThread(new Runnable()   {
	        public void run()  {
	          	AlertDialog alertDialog = new AlertDialog.Builder(activity).create();
	          	alertDialog.setTitle("Error");
	          	alertDialog.setMessage(msg);
	          	alertDialog.setButton(DialogInterface.BUTTON_POSITIVE, "OK",
	          		new DialogInterface.OnClickListener() {
	          			public void onClick(DialogInterface dialog, int i) {}
	          	});
	          	alertDialog.show();
	        }
	    });    
    }

    private void loadInitialSeries() {
	
//		Uri uri = Uri.parse("file:///storage/sdcard0/AudioTime/4007.wav");
		
//    	Uri uri = getIntent().getData();
    	
		plot.setTitle(getString(R.string.app_name) + " " + uri.toString());

		byte data[] = new byte[2];

		WaveHeader waveHeader = new WaveHeader();
		InputStream in = null;
		
		try {
			ContentResolver cr = getContentResolver();
			in = cr.openInputStream(uri);
			waveHeader.read(in);
		} catch (Exception e) {
			errorDialog("WAV file incorrect formatting or not found: " + uri.toString());
			e.printStackTrace();
			return;
		}

		int samplesAvailable = waveHeader.getNumBytes() / 2 ;
		
		if(samplesAvailable > MAXSAMPLES)
			samples = (int)(samplesAvailable/MAXSAMPLES)*MAXSAMPLES;								// Multiple of MAXSAMPLES
		else
			samples =  (int) Math.pow(2.0, (int)(Math.log10(samplesAvailable) / Math.log10(2))); 	// power of 2

		sampleRate = waveHeader.getSampleRate();
		
		int total = 0, n = 0, index = 0;
		int N = Math.min(samples, MAXSAMPLES);
		int y;
		double audioData[] = new double[N];	
		int halfSpectrum = sampleRate / 2;
		double average[] = new double[Math.min(halfSpectrum, N/2)];
		double mag[];
		FFTresult fftResult = null;
	//	Log.d("******loadInitialSeries ",samples+" " +sampleRate+" " + waveHeader.getNumBytes()+" " +N);

		try {
			while (in.read(data, 0, 2) != -1 && total < samples) {
				y = (data[0] & 0xFF) | ((data[1] & 0xFF) << 8); 			// Little endian
				y = y <= 32767 ? y : y - 65535;
				audioData[index] = y;
				if( index == N - 1) {
					mag = (new FFT()).fft(audioData);
					fftResult = new FFTresult(mag, sampleRate, N);
					for(int i=0; i < fftResult.magnitudeFFT.length; i++)
						average[i] = average[i]+fftResult.magnitudeFFT[i];	// Sum normalize spectrum magnitudes
					n++;
					index = 0;
				}
				else	
					index++;
				total++;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		if( fftResult == null ) {
			errorDialog("WAV file incorrect sample number: " + uri.toString());
			return;
		}

		for(int i=0; i < fftResult.magnitudeFFT.length; i++)					// Ensemble spectral averaging
			average[i] = average[i]/n;
		
		fftResult.setMagnitudeFFT(average);
		
		graphSpectrum(fftResult.magnitudeFFT, fftResult.frequencyFFT, sampleRate, samples);
	}
	
	private void graphSpectrum(double []newMagnitudeFFT, int []newFrequencyFFT, int newSampleRate, int newSamples) {

		sampleRate = newSampleRate;
		samples = newSamples;
		frequencyFFT = newFrequencyFFT;
		magnitudeFFT = newMagnitudeFFT;

		if(samples < 16) return;
		
		int spectrumSize = magnitudeFFT.length;
				
		int maxI = 0, minI=0;
		
		for (int i = 0; i < spectrumSize; i++) {
			if (magnitudeFFT[i] > magnitudeFFT[maxI]) 
				maxI = i;
			if (magnitudeFFT[i] < magnitudeFFT[minI])
				minI = i;
			if(frequencyFFT[i] == 0)
				frequencyFFT[i] = 1;												// Remove any 0 frequency
		}				
//		Log.d("******graphSpectrum ",frequencyFFT.length + " " + magnitudeFFT.length);
		
		samplesToDisplay = Math.min(spectrumSize, SPECTRUMSIZE);

		int spectrumStep = spectrumSize/samplesToDisplay;
		
		seriesData = new SimpleXYSeries("Hz");

		seriesData.addLast(frequencyFFT[0], magnitudeFFT[0]);
		
		for (int i = 0; i < spectrumSize; i=i+spectrumStep) {						// Construct spectrum of highest magnitudes within each step size
			int maxJ=i;
			for(int j=i; j<i+spectrumStep && j < spectrumSize; j++) {
				if(magnitudeFFT[j]>magnitudeFFT[maxJ])
					maxJ=j;
			}
			seriesData.addLast(frequencyFFT[maxJ], magnitudeFFT[maxJ]);
		}
		
		/*
		 * freq = i * Fs / N
		 * 
		 * freq = frequency in Hz i = index of peak 
		 * Fs = sample rate (e.g. 44100 Hz or whatever you are using) 
		 * N = size of FFT (e.g. 1024 in your case)
		 */
		plot.setRangeBoundaries(magnitudeFFT[minI], magnitudeFFT[maxI], BoundaryMode.FIXED);
		plot.setDomainBoundaries(seriesData.getX(0).intValue(), seriesData.getX(seriesData.size()-1).intValue(), BoundaryMode.FIXED);
		plot.addSeries(seriesData, seriesFormat);
		
		maxF = frequencyFFT[frequencyFFT.length-1];
		maxX = maxF;
		minF = frequencyFFT[0];

		xMarker = markerAtX(frequencyFFT[maxI], Color.BLACK);
		
//		reloadSeries();
		
//		redrawPlot();
		
//	Log.d("****initialSpectrum",magnitudeFFT.length+" "+samples+" "+sampleRate+" "+seriesData.size()+" "+seriesData.getX(0).intValue() +" "+ seriesData.getX(seriesData.size()-1).intValue()+" "+spectrumSize+" "+spectrumStep);
//		for( int i = 0; i < seriesData.size(); i++)
//			Log.d("\t", seriesData.getX(i)+ " "+seriesData.getY(i));

	}
	
	private void setButtonHandlers() {
		((ImageButton) findViewById(R.id.btnLinearLog))
		.setOnClickListener(new View.OnClickListener() {
			public void onClick(View arg0) {
				if ("isChecked".equals(arg0.getTag())) {
					arg0.setTag("isNotChecked");
					scaleLinear=false;
					((ImageButton) findViewById(R.id.btnLinearLog))
							.setImageResource(R.drawable.log);
					onScaleCheckBoxClicked(false);
				} else {
					arg0.setTag("isChecked");
					scaleLinear=true;
					((ImageButton) findViewById(R.id.btnLinearLog))
							.setImageResource(R.drawable.linear);
					onScaleCheckBoxClicked(true);
				}
			}
		});

		((ImageButton) findViewById(R.id.btnReset))
				.setOnClickListener(new View.OnClickListener() {
					@Override
					public void onClick(View arg0) {
						minF = 1f;
						maxF = maxX;
						reloadSeries();
						redrawPlot();
					}
				});

		((ImageButton) findViewById(R.id.btnLeft))
		.setOnClickListener(new View.OnClickListener() {
			public void onClick(View arg0) {
				maxF=(maxF-minF);
				minF=1f;
				reloadSeries();
				redrawPlot();
			}
		});

		((ImageButton) findViewById(R.id.btnRight))
		.setOnClickListener(new View.OnClickListener() {
			public void onClick(View arg0) {
				minF=maxX-(maxF-minF);
				maxF=maxX;
				reloadSeries();
				redrawPlot();
			}
		});

		((ImageButton) findViewById(R.id.btnHelp))
				.setOnClickListener(new View.OnClickListener() {
					public void onClick(View arg0) {
						AlertDialog.Builder builder = new AlertDialog.Builder(activity);

						final ListView modeList = new ListView(activity);

						ArrayAdapter<String> modeAdapter = new ArrayAdapter<String>(
								activity, android.R.layout.simple_list_item_1,
								android.R.id.text1, DIALOGMENU);
						modeList.setAdapter(modeAdapter);

						builder.setView(modeList);
						final Dialog dialog = builder.create();
						dialog.setCanceledOnTouchOutside(true);

						dialog.show();
						modeList.setOnItemClickListener(new android.widget.AdapterView.OnItemClickListener() {
							@Override
							public void onItemClick(AdapterView<?> parent,
									View view, int position, long id) {
								switch (position) {
								case ABOUT:
									AboutDialog about = new AboutDialog(activity);
									about.setTitle("About AudioSpectrum");
									about.show();
									break;
								case INFO:
									InfoDialog info = new InfoDialog(activity);
									info.setTitle("Info AudioSpectrum");
									info.show();
									break;
								}
								dialog.dismiss();
							}
						});
					}
				});
	}

	private void startReloadThread(){
		threadReload = new Thread(new Runnable() {

		@Override
		public void run() {
				reload();
			}
		}, "AudioSpectrum Reload Thread");

		threadReload.start();
	}
	
	private void reloadSeries() {
					
		long currentTime=System.currentTimeMillis();
		if(currentTime-lastTime < ELAPSEDTIME) return;							// Only reload when visually necessary
		
		lastTime=currentTime;

		if(maxF-minF <= 1.0) return;

		isReloading = true;
		synchronized(plot) {
			plot.notifyAll();
		}
	}
	
	private void reload() {
		int maxI, minI;
		
		while(true)
			synchronized(plot){ 
				while(!isReloading)	
					try {
						plot.wait();
					} catch (InterruptedException e) {}
		
				plot.removeSeries(seriesData);
		
				seriesData = new SimpleXYSeries("Hz");
				
				if(samples < sampleRate) {
					maxI = (int) (frequencyFFT.length * (maxF/sampleRate))*2;
					minI = (int) (frequencyFFT.length * (minF/sampleRate))*2;
				}
				else {
					maxI = (int) maxF;
					minI = (int) minF;					
				}
					
				int spectrumStep = (maxI-minI)/samplesToDisplay;
				if(spectrumStep <= 0)
					spectrumStep=1;

//				Log.d("0 *******reload ", maxF+" "+maxI+" "+minF+" "+minI+" " + (maxI-minI) + " " +frequencyFFT.length+" " +sampleRate+" "+samples+" "+scaleLinear);

				if(scaleLinear)
					seriesData.addLast(frequencyFFT[minI], magnitudeFFT[minI]);
				else 
					seriesData.addLast(Math.log10(frequencyFFT[minI]), magnitudeFFT[minI]);
				
				for (int i = minI; i < maxI; i=i+spectrumStep) {						// Construct spectrum of highest magnitudes within each step size
					int maxJ=i;
					for(int j=i; j<i+spectrumStep && j < maxI; j++)
						if(magnitudeFFT[j]>magnitudeFFT[maxJ])
							maxJ=j;
					if(scaleLinear || frequencyFFT[maxJ] == 0)
						seriesData.addLast(frequencyFFT[maxJ], magnitudeFFT[maxJ]);
					else
						seriesData.addLast(Math.log10(frequencyFFT[maxJ]), magnitudeFFT[maxJ]);
				}
		/*
		 * freq = i * Fs / N
		 * 
		 * freq = frequency in Hz i = index of peak 
		 * Fs = sample rate (e.g. 44100 Hz or whatever you are using) 
		 * N = size of FFT (e.g. 1024 in your case)
		 */

			plot.setDomainBoundaries(seriesData.getX(0), seriesData.getX(seriesData.size()-1), BoundaryMode.FIXED);
//				Log.d("1******setDomainBoundaries start end scaleLinear ",seriesData.getX(0) + " " + seriesData.getX(seriesData.size()-1)+" "+scaleLinear);
//				Log.d("1 ******reload maxF minF maxX ",maxF + " " + minF+" "+maxX);
			plot.addSeries(seriesData, seriesFormat);
			correctMarker();
			isReloading = false;
		}
	}

	private void startRedrawThread(){
		threadRedraw = new Thread(new Runnable() {

		@Override
		public void run() {
				redraw();
			}
		}, "AudioSpectrum Redraw Thread");

		threadRedraw.start();
	}
	
	private void redrawPlot() {
		isDrawing = true;
		synchronized(plot) {
			plot.notifyAll();
		}
	}

	private void redraw() {
		  while(true)
				synchronized(plot){ 
					while(!isDrawing)	
						try {
							plot.wait();
						} catch (InterruptedException e) {}		
					plot.redraw();
					isDrawing = false;
				}
	}

	private void onScaleCheckBoxClicked(boolean checked) {
		
		if(seriesData.size() == 0) return;
		
		int n;
		
		plot.removeMarkers();

		if (checked) {
			n=seriesData.size();
			for (int i = 0; i < n; i++) {
				seriesData.addLast( Math.round(Math.pow(10, seriesData.getX(0).doubleValue())), seriesData.getY(0));
				seriesData.removeFirst();
			}
			graphXLabelFormat.setScaleLinear(true); 
			xMarker=markerAtX(Math.round(Math.pow(10,xMarker.getValue().doubleValue())), Color.BLACK);
		} else {
			n=seriesData.size();
			for (int i = 0; i < n; i++) {
				seriesData.addLast(Math.log10(seriesData.getX(0).doubleValue()+0.1), seriesData.getY(0));
				seriesData.removeFirst();
			}
			graphXLabelFormat.setScaleLinear(false);
			xMarker=markerAtX(Math.log10(xMarker.getValue().doubleValue()+0.1), Color.BLACK);
		}
		plot.setDomainBoundaries(seriesData.getX(0), seriesData.getX(seriesData.size()-1), BoundaryMode.FIXED);
		redrawPlot();
	}

	private void onPlotClicked(PointF point) {
		// make sure the point lies within the graph area. we use gridrect
		// because it accounts for margins and padding as well.
		if (plot.getGraphWidget().getGridRect().contains(point.x, point.y)) {
			Number x = plot.getXVal(point);
			Number y = plot.getYVal(point);

			selection = null;
			double xDistance = 0;
			double yDistance = 0;

			// find the closest value to the selection:
			for (XYSeries series : plot.getSeriesSet()) {
				for (int i = 0; i < series.size(); i++) {
					Number thisX = series.getX(i);
					Number thisY = series.getY(i);
					if (thisX != null && thisY != null) {
						double thisXDistance = LineRegion.measure(x, thisX)
								.doubleValue();
						double thisYDistance = LineRegion.measure(y, thisY)
								.doubleValue();
						if (selection == null) {
							selection = new Pair<Integer, XYSeries>(i, series);
							xDistance = thisXDistance;
							yDistance = thisYDistance;
						} else if (thisXDistance < xDistance) {
							selection = new Pair<Integer, XYSeries>(i, series);
							xDistance = thisXDistance;
							yDistance = thisYDistance;
						} else if (thisXDistance == xDistance
								&& thisYDistance < yDistance
								&& thisY.doubleValue() >= y.doubleValue()) {
							selection = new Pair<Integer, XYSeries>(i, series);
							xDistance = thisXDistance;
							yDistance = thisYDistance;
						}
					}
				}
			}

		} else 
			selection = null;				// if the press was outside the graph area, deselect:

		plot.removeXMarkers();

		if (selection != null)
			markerAtX(selection.second.getX(selection.first), Color.BLACK);

		plot.redraw();
	}
	
	// Definition of the touch states
	static final int NONE = 0;
	static final int ONE_FINGER_DRAG = 1;
	static final int ONE_FINGER_ADJUST = 2;
	static final int TWO_FINGERS_DRAG = 3;
	static int mode = NONE;

	static PointF firstFinger;
	static float distBetweenFingers;

	private void zoom(float scale) {
		float domainSpan = maxF - minF;
		float domainMidPoint = maxF - domainSpan / 2.0f;
		float offset = domainSpan * scale / 2.0f;
		float leftBoundary = 0f; 		
		float rightBoundary = maxX; 	

		minF = domainMidPoint - offset;
		maxF = domainMidPoint + offset;

		if (minF < leftBoundary)
			minF = leftBoundary;
		if (maxF > rightBoundary)
			maxF = rightBoundary;
	}

	private void scroll(float pan) {
		float domainSpan = maxF - minF;
		float step = domainSpan / plot.getWidth();
		float offset = pan * step;
		float leftBoundary = 0f; 
		float rightBoundary = maxX; 

		minF = minF + offset;
		maxF = maxF + offset;

		if (minF < leftBoundary) {
			minF = leftBoundary;
			maxF = leftBoundary + domainSpan;
		} else if (maxF > rightBoundary) {
			maxF = rightBoundary;
			minF = rightBoundary - domainSpan;
		}
	}

	private float spacing(MotionEvent event) {
		float x = event.getX(0) - event.getX(1);
		float y = event.getY(0) - event.getY(1);
		return (float) Math.sqrt(x * x + y * y);
	}

	@Override
	public boolean onTouchEvent(MotionEvent event) {

		if (!plot.getGraphWidget().getGridRect().contains(event.getX(), event.getY())) 
			return super.onTouchEvent(event); 								// Ignore stray touches outside plot area			

		detector.onTouchEvent(event);

		switch (event.getAction() & MotionEvent.ACTION_MASK) {
			case MotionEvent.ACTION_DOWN: 
				firstFinger = new PointF(event.getX(), event.getY());
				mode = ONE_FINGER_DRAG;
				break;
			case MotionEvent.ACTION_UP:
			case MotionEvent.ACTION_POINTER_UP:
				mode = NONE;
				break;
			case MotionEvent.ACTION_POINTER_DOWN: 								// second finger
				distBetweenFingers = spacing(event);							// the distance check is done to avoid false alarms
				if (distBetweenFingers > 5f) {
					mode = TWO_FINGERS_DRAG;
				}
				break;
			case MotionEvent.ACTION_MOVE:
				if (mode == ONE_FINGER_DRAG) {
					PointF oldFirstFinger = firstFinger;
					firstFinger = new PointF(event.getX(), event.getY());
					if (Math.abs(oldFirstFinger.x - firstFinger.x) > 0) {
						scroll(oldFirstFinger.x - firstFinger.x);
						reloadSeries();
						redrawPlot();
					}
				} else if (mode == TWO_FINGERS_DRAG) {
					float oldDist = distBetweenFingers;
					distBetweenFingers = spacing(event);
					zoom(oldDist / distBetweenFingers);
					reloadSeries();
					redrawPlot();
				}
				break;
		}

		return super.onTouchEvent(event);
	}

	private void correctMarker() {
//		Log.d("1******correctMarker ",maxF+" "+minF+" "+seriesData.getX(0).floatValue()+" "+seriesData.getX(seriesData.size()-1).floatValue()+" "+xMarker.getValue().floatValue());
		if (isMarkerOffGraph(xMarker))
			plot.removeMarker(xMarker);

		if (isMarkerOnGraph(xMarker))
			plot.addMarker(xMarker);
	}

	private boolean isMarkerOffGraph(XValueMarker m) {
		if (m == null)
			return false;
		return m.getValue().floatValue() < seriesData.getX(0).floatValue() || m.getValue().floatValue() > seriesData.getX(seriesData.size()-1).floatValue();
	}

	private boolean isMarkerOnGraph(XValueMarker m) {
		if (m == null)
			return false;
		return m.getValue().floatValue() >= seriesData.getX(0).floatValue() && m.getValue().floatValue() <= seriesData.getX(seriesData.size()-1).floatValue();
	}

	@Override
	public boolean onDown(MotionEvent e) {
		// Log.d("---onDown----", e.toString());
		return false;
	}

	@Override
	public boolean onFling(MotionEvent e1, MotionEvent e2, float velocityX,
			float velocityY) {
//		 Log.d("---onFling---", e1.toString() + e2.toString());
		return false;
	}

	@Override
	public void onLongPress(MotionEvent e) {
		// Log.d("---onLongPress---", e.toString());
	}

	@Override
	public boolean onScroll(MotionEvent e1, MotionEvent e2, float distanceX,
			float distanceY) {
		// Log.d("---onScroll---", e1.toString() + e2.toString());
		return false;
	}

	@Override
	public void onShowPress(MotionEvent e) {
		// Log.d("---onShowPress---", e.toString());
	}

	@Override
	public boolean onSingleTapUp(MotionEvent e) {
//		 Log.d("---onSingleTapUp---", e.toString());
		return false;
	}

	@Override
	public boolean onDoubleTap(MotionEvent e) {
		// Log.d("---onDoubleTap---", e.toString());
		return false;
	}

	@Override
	public boolean onDoubleTapEvent(MotionEvent e) {
		// Log.d("---onDoubleTapEvent---", e.toString());
		return false;
	}

	@Override
	public boolean onSingleTapConfirmed(MotionEvent e) {
//		Log.d("---onSingleTapConfirmed", e.toString());
		onPlotClicked(new PointF(e.getX(), e.getY()));
		return false;
	}
	
	private XValueMarker markerAtX(Number xVal, int color) {
		int label = scaleLinear ? xVal.intValue() : (int) Math.round(Math.pow(10, xVal.doubleValue()));
		
		xMarker = new XValueMarker(xVal, 					// xVal to mark
				label + " Hz", new YPositionMetric( 		// object instance to set text positioning on the marker
						PixelUtils.dpToPix(5), 				// 5dp offset
						YLayoutStyle.ABSOLUTE_FROM_TOP), 	// offset origin
				color, 										// line paint color
				color); 									// text paint color

		xMarker.getTextPaint().setTextSize(PixelUtils.dpToPix(14));

		DashPathEffect dpe = new DashPathEffect(new float[] {
				PixelUtils.dpToPix(2), PixelUtils.dpToPix(2) }, 0);

		xMarker.getLinePaint().setPathEffect(dpe);
		xMarker.setTextOrientation(ValueMarker.TextOrientation.VERTICAL);

		plot.addMarker(xMarker);
		plot.redraw();
		return xMarker;
	}
}

class GraphXLabelFormat extends Format {
	private static final long serialVersionUID = 1L;
	private boolean linear = true;

	public GraphXLabelFormat(boolean linear) {
		this.linear = linear;
	}

	public void setScaleLinear(boolean linear) {
		this.linear = linear;
	}

	@Override
	public StringBuffer format(Object object, StringBuffer buffer, FieldPosition field) {
		int parsedInt;
		if (linear)
			parsedInt = (int) Math.round(Double.parseDouble(object.toString()));
		else
			parsedInt = (int) Math.round(Math.pow(10, Double.parseDouble(object.toString())));

		buffer.append(parsedInt);
		return buffer;
	}

	@Override
	public Object parseObject(String string, ParsePosition position) { return null;	}
}



