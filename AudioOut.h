//
//  AudioOut.h
//
// Contains functions to write an array as an audio file.
// OpenAL has been used to play audio out from command line

#ifndef AudioOut_hpp
#define AudioOut_hpp

#include <iostream>
#include <math.h>

#ifdef __APPLE__
	#include "TargetConditionals.h"
    #ifdef TARGET_OS_MAC
		#include <OpenAL/al.h>	// macOS Only
		#include <OpenAL/alc.h>	// macOS Only
    #endif
#elif defined _WIN32 || defined _WIN64
	#include <AL/al.h>
	#include <AL/alc.h>
#endif
// alut.h library deprecated in macOS: but you can use freealut.
//
// In Terminal
// >> brew install freealut
//
// Don't forget to add the AL folder in the header search path
//
// /usr/local/Cellar/freealut/1.1.0/include
//
#include <AL/alut.h>
#define BACKEND "alut"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~WRITING .WAV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
// WAV header format: http://soundfile.sapp.org/doc/WaveFormat/

typedef struct waveFormatHeader {
	char     ChunkId[4];
	uint32_t ChunkSize;
	char     Format[4];
	char     Subchunk1ID[4];
	uint32_t Subchunk1Size;
	uint16_t AudioFormat;
	uint16_t NumChannels;
	uint32_t SampleRate;
	uint32_t ByteRate;
	uint16_t BlockAlign;
	uint16_t BitsPerSample;
	char     SubChunk2ID[4];
	uint32_t Subchunk2Size;
} waveFormatHeader_t;


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Define a standard 16bit PCM .wav file header.q

waveFormatHeader_t * basicHeader(){
	waveFormatHeader_t * wh = (waveFormatHeader_t *)malloc(sizeof(waveFormatHeader_t));
	memcpy(wh->ChunkId, &"RIFF", 4);
	memcpy(wh->Format, &"WAVE", 4);
	memcpy(wh->Subchunk1ID, &"fmt ", 4); //notice the space at the end!
	wh -> Subchunk1Size = 16;
	return wh;
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// Fill in the gaps of basicHeader() to create a header for a Stereo .wav file

waveFormatHeader_t * stereo16bitWaveHeader(double SampRate){
	waveFormatHeader_t * wh = basicHeader();
	wh->AudioFormat = 1;
	wh->NumChannels = 2;
	wh -> SampleRate = (uint32_t)SampRate;
	wh -> BitsPerSample = 16;
	wh -> ByteRate = wh->NumChannels * wh -> SampleRate * wh -> BitsPerSample/8;
	wh -> BlockAlign = wh->NumChannels * wh -> BitsPerSample/8;
	memcpy(wh->SubChunk2ID, &"data", 4);
	return wh;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\



// Fill in the gaps of basicHeader() to create a header for a Mono .wav file

waveFormatHeader_t * mono16bitWaveHeader(double SampRate){
	waveFormatHeader_t * wh = basicHeader();
	wh->AudioFormat = 1;
	wh->NumChannels = 1;
	wh -> SampleRate = (uint32_t)SampRate;
	wh -> BitsPerSample = 16;
	wh -> ByteRate = wh->NumChannels * wh -> SampleRate * wh -> BitsPerSample/8;
	wh -> BlockAlign = wh->NumChannels * wh -> BitsPerSample/8;
	memcpy(wh->SubChunk2ID, &"data", 4);
	return wh;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

// numberOfFrames in this case referes to the number of samples in the incoming
// array.

void setLengthForWaveFormatHeader(waveFormatHeader_t * wh, size_t numberOfFrames){
	wh -> Subchunk2Size = (uint32_t)numberOfFrames * wh->NumChannels * wh->BitsPerSample/8;
	wh -> ChunkSize = 36 + wh -> Subchunk2Size;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


waveFormatHeader_t * stereo16bitWaveHeaderForLength(size_t numberOfFrames,double SampRate){
	waveFormatHeader_t * wh = stereo16bitWaveHeader(SampRate);
	setLengthForWaveFormatHeader(wh, numberOfFrames);
	return wh;
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


waveFormatHeader_t * mono16bitWaveHeaderForLength(size_t numberOfFrames,double SampRate){
	waveFormatHeader_t * wh = mono16bitWaveHeader(SampRate);
	setLengthForWaveFormatHeader(wh, numberOfFrames);
	return wh;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\



void normaliseBuffer(double *audio, int NF) {
	
	int n;
	double temp;
	
	// Find max abs sample
	double maxy = 0.0;
	
	for(n=0;n<NF;n++)
	{
		if(fabs(audio[n])>maxy) maxy = fabs(audio[n]);
		
	}
	
	// Normalise
	if(maxy>0.00001){
		for(n=0;n<NF;n++)
		{
			temp     = audio[n];
			audio[n] = temp/maxy;
		}
	}
	
	// Smooth last 500 samples
	if(NF>501){
		double inc  = 1.0/500.0;
		double ramp = 1.0;
		for(n=NF-501;n<NF;n++)
		{
			audio[n] *= ramp;
			if(ramp>0) ramp-=inc;
		}
	}
	
	printf("Normalised by : %.5f\n", maxy);
	
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\



size_t writeWaveHeaderToFile(waveFormatHeader_t * wh, FILE * file){
	return fwrite(wh, sizeof(waveFormatHeader_t), 1, file);
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


static void writeWavMS(double *audio,const char outputFile[], int NumberOfSamples, double SampRate){
	
	normaliseBuffer(audio ,NumberOfSamples);
	
	FILE * file;
	file = fopen(outputFile, "w");
	waveFormatHeader_t * wh = mono16bitWaveHeaderForLength(NumberOfSamples,SampRate);
	writeWaveHeaderToFile(wh, file);
	free(wh);
	int16_t sdata;
	double amp = 32000.0;
	
	int i;
	
	for(i=0;i<NumberOfSamples;i++){
		
		sdata = (int16_t)(audio[i]*amp);  //set sdata to PCM 16-bit
		fwrite(&sdata, sizeof(int16_t), 1, file);
		
	}
	fclose(file);
	printf("%d samples written to %s\n", NumberOfSamples,outputFile);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PLAYING .WAV~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\


static void list_audio_devices(const ALCchar *devices)
{
	const ALCchar *device = devices, *next = devices + 1;
	size_t len = 0;
	
	fprintf(stdout, "Devices list:\n");
	fprintf(stdout, "----------\n");
	while (device && *device != '\0' && next && *next != '\0') {
		fprintf(stdout, "%s\n", device);
		len = strlen(device);
		device += (len + 1);
		next += (len + 2);
	}
	fprintf(stdout, "----------\n");
}

#define TEST_ERROR(_msg)		\
error = alGetError();		\
if (error != AL_NO_ERROR) {	\
fprintf(stderr, _msg "\n");	\
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

static inline ALenum to_al_format(short channels, short samples)
{
	bool stereo = (channels > 1);
	
	switch (samples) {
		case 16:
			if (stereo)
				return AL_FORMAT_STEREO16;
			else
				return AL_FORMAT_MONO16;
		case 8:
			if (stereo)
				return AL_FORMAT_STEREO8;
			else
				return AL_FORMAT_MONO8;
		default:
			return -1;
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

#ifdef AL_ALUT_H
// This function reads in a .wav file and plays it out using OpenAL
void playWavMS(char *inputfname){
	
	ALboolean enumeration;
	const ALCchar *defaultDeviceName = NULL;
	ALCdevice *device;
	ALvoid *data;
	ALCcontext *context;
	ALsizei size, freq;
	ALenum format;
	ALuint buffer, source;
	ALCenum error;
	ALint source_state;

	fprintf(stdout, "Using " BACKEND " as audio backend\n\n");

	enumeration = alcIsExtensionPresent(NULL, "ALC_ENUMERATION_EXT");
	if (enumeration == AL_FALSE)
		fprintf(stderr, "enumeration extension not available\n");

	list_audio_devices(alcGetString(NULL, ALC_DEVICE_SPECIFIER));

	if (!defaultDeviceName)
		defaultDeviceName = alcGetString(NULL, ALC_DEFAULT_DEVICE_SPECIFIER);

	device = alcOpenDevice(defaultDeviceName);
	if (!device) {
		
		TEST_ERROR("unable to open default device\n");
	}

	fprintf(stdout, "Device: %s\n", alcGetString(device, ALC_DEVICE_SPECIFIER));

	alGetError();

	context = alcCreateContext(device, NULL);
	if (!alcMakeContextCurrent(context)) {
		fprintf(stderr, "failed to make default context\n");

	}

	TEST_ERROR("make default context");
	
	alGenSources((ALuint)1, &source);
	TEST_ERROR("source generation");
	
	alGenBuffers(1, &buffer);
	TEST_ERROR("buffer generation");
	
	alutLoadWAVFile(inputfname, &format, &data, &size, &freq);
	TEST_ERROR("loading wav file");
	
	alBufferData(buffer, format, data, size, freq);
	TEST_ERROR("buffer copy");
	
	alSourcei(source, AL_BUFFER, buffer);
	TEST_ERROR("buffer binding");
	
	alSourcePlay(source);
	TEST_ERROR("source playing");
	
	alGetSourcei(source, AL_SOURCE_STATE, &source_state);
	TEST_ERROR("source state get");
	while (source_state == AL_PLAYING) {
		alGetSourcei(source, AL_SOURCE_STATE, &source_state);
		TEST_ERROR("source state get");
	}
	
	/* exit context */
	alDeleteSources(1, &source);
	alDeleteBuffers(1, &buffer);
	device = alcGetContextsDevice(context);
	alcMakeContextCurrent(NULL);
	alcDestroyContext(context);
	alcCloseDevice(device);
	
}
#else
void playWavMS(char *inputfname){
	printf("ALUT not defined\n");
}
#endif
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\

void weirdPrint(const char input[]){
	printf("%s\n", input);
}

#endif /* AudioOut_hpp */
