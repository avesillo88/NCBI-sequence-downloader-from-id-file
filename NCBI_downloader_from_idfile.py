#!/usr/bin/python3

#  This script is intended to download amino-acid or nucleotide sequences from a list of identifiers 
#  stored in a file. For this purpose, the NCBI E-utilities Efecth and Epost are used. 

#  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY

#  Usage: NCBI_downloader_from_idfile.py "input_file" "output_file" batch_size "db" "ret_type" "api_key" "email"

import requests
from tqdm import tqdm
import sys
import socket
from os import stat
import time
from Bio import Entrez
import urllib3
import os
import re



class NCBISequenceDownloader:
	def __init__(self, input_file: str, output_file: str, batch_size: int, db: str, ret_type: str, api_key: str, email: str):
		self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
		self.input_file = input_file
		self.batch_size = int(batch_size)
		self.db = db
		self.ret_type = ret_type
		self.output_file = output_file
		self.api_key = api_key
		self.email = email
		self.ids_for_download = self.check_file()
	
	def check_file(self):
		# Check if the output file exists and has a size greater than zero
		if os.path.exists(self.output_file) and os.path.getsize(self.output_file) > 0:
			
			print(f"{self.output_file} exists and has a size greater than zero.")
			
			# Loop until a valid response of "1", "2" or "3" is entered
			while True:
				response = input(f"\nThe output file seems to has some information. If it is not empty, " 
							 " then it may cause errors during the running."
							 " Do you want to remove the file (1), check the output file for protein IDs (2) or exit (3)? (1/2/3) ")
				# Remove the file if the user enters "1"
				if response == "1":
					os.remove(self.output_file)
					ids_for_download = self.get_id_list_to_download()
					break
				# Check output file for IDs if the user enters "2" and compare with the input file
				elif response == "2":
					initial_ids = set(self.get_id_list_to_download())
					ids_already_download = self.get_ids_already_downloaded()
					ids_for_download = list(initial_ids.difference(ids_already_download))
					break
				#Finish the programm without downloading any sequence
				elif response == "3":
					self.finish(False)
					break
				# Print error message if the user enters an invalid response
				else:
					print("Invalid response. Please enter 1, 2 or 3.")
		# If the file doesn't exist or has a zero size, get the list of IDs for download from the input file
		else:
			ids_for_download = self.get_id_list_to_download()
		print(f"\nStarting the download of: {len(ids_for_download)} sequences\n")
		return ids_for_download

			
	
	def get_id_list_to_download(self):
		# Read the file containing the protein IDs, removing any duplicates
		with open(self.input_file) as input_file:
			lines = set(input_file.read().splitlines())
			protein_ids = list(lines)

		return protein_ids
	
	def get_ids_already_downloaded(self):
		with open(self.output_file) as output_file:
			protein_ids = []
			# Iterate over each line in the file
			for line in output_file:
				# Check if the line starts with ">"
				if line.startswith(">"):
					# check if the current line contains the symbol "|" in the first 15 characters
					if "|" in line[:15]:
						# If it contains the "|" symbol, we use this as a delimiter for splitting.
						protein_id = line.split("|")[1]
						if protein_id == '':
							protein_id = line.split("|")[2]
					else:
						# If it does not contain the "|" symbol, we use the space as a delimiter for splitting.
						protein_id = line.split(" ")[0] 

						# Remove the ">" character
						protein_id = protein_id[1:]

					# Add the protein ID to the list
					protein_ids.append(protein_id)
		protein_ids = list(set(protein_ids))
		obtained_ids = set(protein_ids)
		
		
		return obtained_ids

	def download_sequences(self, protein_ids):

		
		error_flag = False

		# Open the file to write the sequences to
		with open(self.output_file, "a") as f:
			# Use tqdm to create a progress bar
			for i in tqdm(range(0, len(protein_ids), self.batch_size)):
				
				# Get the current batch of IDs
				ids = protein_ids[i : i + self.batch_size]

				# Use epost to get the WebEnv and QueryKey for the current batch of IDs
				epost_params = {
					"db": self.db,
					"id": ",".join(ids),
					"api_key": self.api_key,
					"email": self.email,
				}
				
				# Use a try-except block to handle the error
				try:
					epost_response = requests.post(self.base_url + "epost.fcgi", data=epost_params)
					epost_response.raise_for_status()

					# Get the WebEnv and QueryKey, handling different error if they occur
					try:
						webenv = epost_response.text.split("<WebEnv>")[1].split("</WebEnv>")[0]
						query_key = epost_response.text.split("<QueryKey>")[1].split("</QueryKey>")[0]


					except IndexError as error:
						print(f"An error occurred when parsing the response from epost:", error, ". Moving to next set of protein IDs.\n")
						error_flag = True
						# Skip the rest of the current iteration and go to the next one
						continue
					except requests.exceptions.MaxRetryError as e:
						print("Received max retries exceeded with url error. Trying again...")
						time.sleep(10)
						continue
					except requests.exceptions.ConnectionError as e:
						print("Failed to connect:", e, "Trying again...")
						time.sleep(10)
						continue
					except urllib3.exceptions.NewConnectionError:
						print("Failed to connect:", e, "Trying again...")
						time.sleep(10)
						continue
					except: 
						print("Got an empty reply or a Timeout from NCBI."
							  " Let's wait 8'' and try again.")
						#attempt += 1
						sleep(8)

				except socket.gaierror as error:
					print("An error occurred when calling epost:", error)
					error_flag = True
					time.sleep(2)
					# Skip the rest of the current iteration and go to the next one
					continue
					
				except requests.exceptions.HTTPError as err:
					# If the error is a 400 Bad Request, print a message and try again
					#if err.response.status_code == 400:
					print("Received HTTP 400 Bad Request error. Trying again...")
					#attempts += 1
					time.sleep(2)
					error_flag = True
					#continue
					
				except requests.exceptions.ConnectionError as e:
					print("Failed to connect:", e, "Trying again...")
					time.sleep(10)
					continue
				
				except: 
					print("Got an empty reply or a Timeout from NCBI."
						  " Let's wait 8'' and try again.")
					time.sleep(8)
					continue


				# Use efetch to download the sequences
				efetch_params = {
					"db": self.db,
					"query_key": query_key,
					"WebEnv": webenv,
					"rettype": self.ret_type,
					"api_key": self.api_key,
					"email": self.email,
					}
				# Use a try-except block to handle the error
				try:
					efetch_response = requests.post(self.base_url + "efetch.fcgi", data=efetch_params)
					efetch_response.raise_for_status()
				except socket.gaierror as error:
					print("An error occurred when calling efetch:", error)
					error_flag = True
					self.updating_ids_for_download([], 0)
					continue
				except requests.exceptions.HTTPError as error:
					print("An error occurred when calling efetch:", error)
					error_flag = True
					self.updating_ids_for_download([], 0)
					continue

				# Write the sequences to the file
				f.write(efetch_response.text)

				# Sleep for a short time to avoid overloading the NCBI servers
				time.sleep(1)


			
		obtained_ids = self.get_ids_already_downloaded()
		
		return obtained_ids, error_flag
	
	def updating_ids_for_download(self, obtained_ids):


		# Get a list with the total IDs to download
		all_ids_to_download = set(self.ids_for_download)

		
		# Create a set containing the IDs from both files
		ids_left_to_download = list(all_ids_to_download.difference(obtained_ids))
		num_ids_left_to_download = len(ids_left_to_download)

		if num_ids_left_to_download != 0:
			if self.db == 'protein':
				# Dealing with PDB ids
				protein_ids_copy = ids_left_to_download.copy()
				# Open the FASTA file
				with open(self.output_file, 'r') as f:
					# Read line by line
					for line in f:
						# If the line starts with '>'
						if line.startswith('>'):
							# Iterate through the protein id list
							for id_ in protein_ids_copy:
								# Search for the protein identifier in the description line
								match = re.search(r'pdb\|' + id_[:4] + '\|', line)
								if match:
									# Remove the protein identifier from the list
									ids_left_to_download.remove(id_)
				num_ids_left_to_download = len(ids_left_to_download)					
			print ("After this round, the following number of sequences remain to be downloaded: ", num_ids_left_to_download, "sequences")
			
			return False, ids_left_to_download
		else:		 
			
			# finish with a message
			return True, ids_left_to_download


	def efetch_downloader(self, ncbi_ids):
		# Set the email address that will be used for Entrez queries
		Entrez.email = self.email
		Entrez.api_key = self.api_key

		# Set the number of sequences to download in each batch
		batch_size = 10

		# Download the sequences in batches of 9 using the `efetch` function
		with open(self.output_file, "a") as file:
			for i in tqdm(range(0, len(ncbi_ids), batch_size)):
				# Get the next batch of IDs
				batch_ids = ncbi_ids[i:i+batch_size]
				# Download the sequences using the `efetch` function
				seq_download_result = Entrez.efetch(
					db=self.db, rettype=self.ret_type, retmode="text",
					id=",".join(batch_ids), retmax=batch_size
				)
				file.write(seq_download_result.read())
				time.sleep(2)
				
		
		
		protein_ids  = self.get_ids_already_downloaded()
		
		return protein_ids
			
	def main_organizer(self, error_flag = False):
		

		protein_ids = self.ids_for_download
		
		print(f"Starting the first approach to download all sequences with Epost")
		while not error_flag:
			
			protein_ids, error_flag = self.download_sequences(protein_ids)
			success, protein_ids = self.updating_ids_for_download(protein_ids)

			if success:
				self.finish(True, "")
				#break
		else:
			print(f"\nStarting the second approach to download all sequences with Efetch. "
				"Trying to download ", len(protein_ids), "sequences.")


			protein_ids = self.efetch_downloader(protein_ids)
			second_attemp, protein_ids = self.updating_ids_for_download(protein_ids)
			
			if second_attemp: 
				self.finish(True, "")
			else:
				print(f"\nProgram finished without downloading all sequences."
					" It seems some accession numbers cannot be download."
					" Chech manually if they still exist. The ones that"
					" couldn't been downloaded were: ")
				print(protein_ids)
				
				self.finish(False)
		
		
	def finish(self, success, msg=""):
		if success is True:
			print("All sequences were downloaded correctly. Good!")
			sys.exit("Program finished without error.")

		else:
			sys.exit("Program finished with some failures.\n" + msg)

def main():
	dler = NCBISequenceDownloader(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
	dler.main_organizer()


if __name__ == '__main__':
	main()

