/* This program is used to control a Fiber Photometry (FIP) setup. It is part of FIPster....
 *  ....
 *  ...
 *  ...
 */


// Hardware layout
int Camera = 10;
int LED470 = 51;
int LED405 = 50;
int alarm_LED = 13;

// Session-specific variables (default values)
int nr_channels = 2;    //# number of channels
int framerate = 20;     //Hz camera framrate (so divide by 2 if 2 channels)
int tech_delay = 1;     //ms delay caused by Arduino (very small)

// Global variables used for the session (no need to be global....)
unsigned long c_time;   //ms the current estimated time (from counting the loops)
unsigned long s_time;   //ms the time acquisition started
int channels = 2;       //# of channels
int c_delay;            //ms the current delay
long t_delay;           //us tracking average delay
int interval;           //ms interval for 1 channel (includes tech delay)
int o_interval;         //ms interval for 1 loop (calculated from framerate, rounded to 1ms)
int v_interval;         //us variable interval used for time correction
unsigned long loop_counter;       //# counts the loops (used for time correction)

// Global variables used for serial communication
int noun_byte;          //Subject of the instruction
int verb_byte;          //The instruction

// Other variables
int result;             //Check outcome of functino calls


int acquisition(){
  // Controls the actual acquisition

  // Calculating correct interval
  o_interval = 1000/framerate * channels;
  interval = (o_interval/channels - tech_delay) - 5; // in ms (5ms is the camera pulse)

  /* This is not super ellegant, but this is the delay of a the Matlab
   *  function that is used to start the DAQ. To be sure: it's the delay
   *  of a prepared DAQ session on the daq USB-6001. It might be
   *  different on a different DAQ. We incorporate this delay here to
   *  make sure the DAQ AI acquisition and the Arduino are perfectly
   *  alligned.
   */
  delay(75);

  // Start variables
  loop_counter = 0;
  s_time = millis();
  v_interval = 750;
  bool bool_running = true;
  

  while(bool_running){

    // Check if Matlab has something to say
     if(Serial.available() > 1){

      noun_byte = Serial.read();
      verb_byte = Serial.read();

      switch(noun_byte){
        
        // noun 9 is acquisition
        case 9: if(verb_byte==1){
                  // verb 1 means stop
                  bool_running = false;
                  Serial.println("Stopping acquisition on Arduino");
                }else if(verb_byte==3){
                  //verb 3 means request info
                  Serial.println("Will print data here in the future");
                }
                break;
                
        // Other input (not recognized)        
        default: break;
      }
     }

    // Allign the loop with real time
    c_time = loop_counter *  o_interval + s_time;
    while(millis()<c_time){
      //delayMicroseconds(1);
      digitalWrite(alarm_LED, HIGH);
    }
    digitalWrite(alarm_LED, LOW);

    // The actual loop
    digitalWrite(LED405, HIGH);
    digitalWrite(Camera, HIGH);
    delayMicroseconds(5000);
    digitalWrite(Camera, LOW);
    delay(interval);
    digitalWrite(LED405, LOW);
    digitalWrite(LED470, HIGH);
    digitalWrite(Camera, HIGH);
    delayMicroseconds(5000);
    digitalWrite(Camera, LOW);
    delay(interval);
    digitalWrite(LED470, LOW);

    // Error handeling if tech delay is to short and Arduino falling behind
    // TO DO TO DO TO DO TO DO TO DO TO DO
  
    // Increment loop counter
    loop_counter++;
  }
  return 1;
}

void setup() {
  // Setting up the hardware
  pinMode(Camera, OUTPUT);
  pinMode(LED470, OUTPUT);
  pinMode(LED405, OUTPUT);
  pinMode(alarm_LED, OUTPUT);

  // Making Serial connection
  Serial.begin(9600);
}

void loop() {
  
  // Check if there is serial input from the computer
  if(Serial.available() > 1){

    /* Ther are probably more standard ways to mediate serial communication,
     *  but hey, I got this from a youtube video on the board computer of the
     *  Saturn rockets that took people to the moon, so it seemed cool to try
     *  it out.
     */
    noun_byte = Serial.read();
    verb_byte = Serial.read();

    switch(noun_byte){

      // Just want to see if we are connected
      case 1: Serial.println("Arduino is connected and running.");
              break;

      // Set the framerate        
      case 2: framerate = verb_byte; //should do some error handeling....
              Serial.print("Framerate set to: ");
              Serial.println(verb_byte);
              break;

      // Set the tech delay        
      case 3: tech_delay = verb_byte; //
              Serial.print("Tech delay set to: ");
              Serial.println(verb_byte);
              break;
              
      // Start acquisition        
      case 9: Serial.println("Start acquisition on Arduino.");
              result = acquisition();
              break;

      // Unrecognized command        
      default:  Serial.print("Unrecognized command: ");
                Serial.print(noun_byte);
                Serial.print("   ");
                Serial.println(verb_byte);
                break;
    }
    
  }

}
