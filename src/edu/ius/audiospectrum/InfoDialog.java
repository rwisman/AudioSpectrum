package edu.ius.audiospectrum;

import android.app.Dialog;
import android.content.Context;
import android.graphics.Color;
import android.os.Bundle;
import android.widget.TextView;

public class InfoDialog extends Dialog{

	public InfoDialog(Context context) {
		super(context);
	}
	
	/**
     * This is the standard Android on create method that gets called when the activity initialized.
     */
	@Override
	public void onCreate(Bundle savedInstanceState) {
		setCanceledOnTouchOutside(true);
		setContentView(R.layout.info);
		TextView tv = (TextView)findViewById(R.id.info_text);
		tv.setBackgroundColor(Color.WHITE);		
		tv.setTextColor(Color.BLACK);		
	}
}
